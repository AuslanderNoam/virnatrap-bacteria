from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Conv1D, MaxPooling1D, Flatten, BatchNormalization, Embedding
from tensorflow.keras.metrics import AUC, CategoricalAccuracy
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import Nadam, SGD
from tensorflow.keras.callbacks import Callback #, ModelCheckpoint, EarlyStopping, 
from tensorflow.keras.utils import to_categorical
#from tensorflow.keras.models import load_model
#from sklearn.metrics import confusion_matrix, classification_report
#import tensorflow.keras.backend as K
#import tensorflow as tf

import numpy as np
import argparse, sys, os, random
#from Bio import SeqIO
#from tqdm import tqdm
os.system("nvidia-smi")


def get_model(input_shape, hp):
    #print(hp)
    model = Sequential()
    # Convolutional Layers
    #Doesn't allow tuning: different #s of filters/layer, different dropout
    #INIT = "he_normal"
    model.add(Embedding(4, 25, input_length=input_shape[0]))
    widths = [int(x) for x in hp["filtwidth"].split(",")]
    while len(widths) < hp["convlayers"]:
        widths.append(widths[-1])
    for i in range(hp["convlayers"]): #Number of convolutional layers
        model.add(Conv1D(filters = hp["filters"], kernel_size = widths[i],
                         activation='relu', kernel_regularizer = l2(l=hp["l2norm"]),
                         #kernel_initializer=INIT, bias_initializer=INIT, 
                         #input_shape=input_shape, 
                         padding="same"))
        model.add(Dropout(rate = hp["dropout"]))
    if hp["convlayers"] > 0:
        # Pool
        model.add(MaxPooling1D(pool_size=hp["pool"], strides=1, padding="same"))
        #Batch Norm
        model.add(BatchNormalization(scale=False))
    #Flatten
    model.add(Flatten())   
    # Linear Layers
    linunits = [int(x) for x in hp["linunits"].split(",")]
    while len(linunits) < hp["linlayers"]:
        linunits.append(linunits[-1])
    for i in range(hp["linlayers"]):
        model.add(Dense(units = linunits[i], activation = 'relu',
                        #kernel_initializer=INIT, bias_initializer=INIT,
                        kernel_regularizer = l2(l=hp["l2norm"])))
        model.add(Dropout(rate = hp["dropout"]))
    # output layer
    model.add(Dense(units = 3, activation = 'softmax',
                  #kernel_initializer=INIT, bias_initializer=INIT,
                  kernel_regularizer = l2(l=hp["l2norm"])))
    myoptimizer = Nadam(learning_rate=hp["lr"], beta_1=hp["beta"]) #(default beta 0.9)
    model.compile(optimizer=myoptimizer,
                  loss="categorical_crossentropy",
                  metrics=[CategoricalAccuracy(name="accuracy"),
                          AUC(name='auprc_0', curve='PR', multi_label=True, label_weights = [1, 0, 0]),
                          AUC(name='auprc_1', curve='PR', multi_label=True, label_weights = [0, 1, 0]),
                          AUC(name='auprc_2', curve='PR', multi_label=True, label_weights = [0, 0, 1])]
                  )
    #Those AUC metrics give the 1 vs. rest AUPRC for each class of examples. It's a bit inefficient since
    #it basically computes AUPRC for all three classes three times and takes a weighted average w/different weights.
    model.summary()
    return model

#This callback monitors the minimum of CLASS APURCs for early stopping.
#This requires that the minimum class AUPRC improve. An alternative would be to require that
#ALL class AUPRCs improv, but that feels a bit too harsh.
class CallbackClassAUPRC(Callback):
    def __init__(self, patience=0, restore_best_weights = True, save_best_model = False, early_stopping = False, min_delta = 0, model_name = ""):
        super(CallbackClassAUPRC, self).__init__()
        self.patience = patience
        self.restore_best_weights = restore_best_weights
        self.best_weights = None
        self.min_delta = min_delta
        self.early_stopping = early_stopping
        self.save_best_model = save_best_model
        self.model_name = model_name
        
    def on_train_begin(self, logs=None):
        # The number of epoch it has waited when loss is no longer minimum.
        self.wait = 0
        # The epoch the training stops at.
        self.stopped_epoch = 0
        self.best_epoch = 0
        # Initialize the best as 0.
        self.best_min_auprc = 0

    def on_epoch_end(self, epoch, logs=None): 
        min_auprc = min(logs.get('val_auprc_0'), logs.get('val_auprc_1'), logs.get('val_auprc_2'))
        if min_auprc > self.best_min_auprc + self.min_delta:
            self.best_min_auprc = min_auprc
            self.best_epoch = epoch
            self.wait = 0
            if self.restore_best_weights:
                # Record the best weights if current results is better.
                self.best_weights = self.model.get_weights()
            if self.save_best_model:
                self.model.save(self.model_name)
        else:
            self.wait += 1
            if self.early_stopping and self.wait >= self.patience:
                self.stopped_epoch = epoch
                self.model.stop_training = True

    def on_train_end(self, logs=None):
        if self.stopped_epoch > 0:
            print("Epoch %05d: early stopping" % (self.stopped_epoch + 1))
        if self.restore_best_weights:
            print(f'Restoring model weights from the end of the best epoch: {self.best_epoch + 1}.')
            self.model.set_weights(self.best_weights)
        else:
            print(f'The best epoch was: {self.best_epoch + 1}')

def train_model(x_train, y_train, x_valid, y_valid, hp, model):
    # For now, no class weighting
    #total = y_train.shape[0]
    #weight_for_0 = (1 / np.sum(y_train==0))*(total)/2.0
    #weight_for_1 = (1 / np.sum(y_train==1))*(total)/2.0
    #class_weight = {0: weight_for_0, 1: weight_for_1}
    # An epoch is calculated by dividing the number of training images by the batchsize
    
    #y_train_int = np.argmax(y_train, axis=-1)
    #y_train_int[y_train_int == 2] = 1
    #y_train_int[y_train_int == 0] = 2

    #For now, no early saving. If desired, should probably add to custom callback.
    #model_checkpoint_callback = ModelCheckpoint(filepath=hp["name"]+".h5", save_weights_only=False,
    #                                            monitor='val_loss', mode='min', save_best_only=True)
    #early_stop_callback = EarlyStopping(monitor="val_loss", mode="min", patience=5, verbose=1)
    hist = model.fit(x_train, y_train, batch_size = hp["batch"], epochs = hp["epoch"],
                     verbose = 2, #sample_weight = y_train_int,
                     validation_data=(x_valid, y_valid), 
                     callbacks = [CallbackClassAUPRC(restore_best_weights = False, save_best_model = True, model_name = hp["name"])])
    return hist


def load_data(f):
    #Convert from one-hot
    foo = np.load(f)
    x = np.argmax(foo["x"], axis=-1)
    #x = foo["x"]
    #Convert int labels to 1-hot
    y = to_categorical(foo["y"])
    return x, y
    

def main(hp):
    model_name = hp["name"] if hp["name"] else "tmp.h5"
    # load data
    x_train, y_train = load_data(hp["train"])
    if hp["sample"]:
        n = x_train.shape[0]
        indices = random.sample(range(n), k = n // 20)
        x_train = x_train[indices]
        y_train = y_train[indices]
    x_valid, y_valid =  load_data(hp["valid"])
    
    #K.clear_session()  
    # train model
    model = get_model(x_train.shape[1:], hp)
    train_model(x_train, y_train, x_valid, y_valid, hp, model)
    #model.save(model_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--train", type=str, help="Preprocessed Training Examples")
    parser.add_argument("-v", "--valid", type=str, help="Preprocessed Validation Examples")
    parser.add_argument("-b", "--batch", type=int, help="Batch size", default=4096)
    parser.add_argument("-e", "--epoch", type=int, help="Number of epochs")
    parser.add_argument("-r", "--lr", type=float, help="Learning rate")
    parser.add_argument("-2", "--l2norm", type=float, help="L2 normalization", default=0.01)
    parser.add_argument("-f", "--filters", type=int, help="Filters / convolutional layer", default=64)
    parser.add_argument("-u", "--linunits", type=str, help="Units / fully connected layer; comma-separated for individual-layer control", default="64")
    parser.add_argument("-l", "--linlayers", type=int, help="# fully connected layers", default=1)
    parser.add_argument("-p", "--dropout", type=float, help="Dropout", default=0.1)
    parser.add_argument("-o", "--pool", type=int, help="Width of pool", default=5)
    parser.add_argument("-w", "--filtwidth", type=str, help="Width of convolutional filters")
    #parser.add_argument("-m", "--model", type=str, help="Pre-trained TensorFlow model")
    parser.add_argument("-n", "--name", type=str, help="Output model name")
    parser.add_argument("-c", "--convlayers", type=int, help="# convolutional layers", default=1)
    parser.add_argument("-m", "--beta", type=float, help="Momentum decay", default=0.9)
    parser.add_argument("--sample", action="store_true", help="Sample 5% of input?")

    options, args = parser.parse_known_args()
    if (len(sys.argv)==1):
        parser.print_help(sys.stderr)
        sys.exit(1)
    main(vars(options))
