
import sys
import numpy as np
from tensorflow.keras.models import load_model
from sklearn.metrics import confusion_matrix
from tensorflow.keras.metrics import AUC
#sys.path.append("/wistar/auslander/Daniel/scripts")
from cnn3class import load_data
import random

def metrics(labels, preds):
    m = AUC(curve="PR", multi_label=True, label_weights = [1, 0, 0])
    m.update_state(labels, preds)
    print("Virus AUPRC:", m.result().numpy())
    m = AUC(curve="PR", multi_label=True, label_weights = [0, 1, 0])
    m.update_state(labels, preds)
    print("Human AUPRC:", m.result().numpy())
    m = AUC(curve="PR", multi_label=True, label_weights = [0, 0, 1])
    m.update_state(labels, preds)
    print("Bacteria AUPRC:", m.result().numpy())

    print("confusion matrix:")
    print(confusion_matrix(labels.argmax(axis=1), preds.argmax(axis=1)))


valid_x, valid_y = load_data(sys.argv[2])
model = load_model(sys.argv[1])

#print("Original performance:")
#ps = model.predict(valid_x)
#metrics(valid_y, ps)

mut_x = valid_x.copy()
for k in range(1,4):
    for i in range(mut_x.shape[0]):
        j = random.randrange(mut_x.shape[1])
        while mut_x[i][j] != valid_x[i][j]:
            j = random.randrange(mut_x.shape[1])
        mut_x[i][j] = (mut_x[i][j] + random.randrange(1,4)) % 4
    print("Performance with", k, "mutations:")
    ps = model.predict(mut_x)
    metrics(valid_y, ps)
            
    

