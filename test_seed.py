
import sys
import numpy as np
from tensorflow.keras.models import load_model
from sklearn.metrics import confusion_matrix, precision_recall_curve
#sys.path.append("/wistar/auslander/Daniel/scripts")
from cnn3class import load_data
import matplotlib.pyplot as plt


test_x, test_y = load_data(sys.argv[2])
y_labels = test_y.argmax(axis=1)

#vir_idx = np.where(y_labels == 0)[0]
#hum_idx = np.where(y_labels == 1)[0]
#bac_idx = np.where(y_labels == 2)[0]

#rng = np.random.default_rng(0)
#k=1000
#vir_idx_sample = rng.choice(vir_idx, size=k, replace=False)
#hum_idx_sample = rng.choice(hum_idx, size=k, replace=False)
#bac_idx_sample = rng.choice(bac_idx, size=k, replace=False)#

#idx_sample = np.concatenate((vir_idx_sample, hum_idx_sample, bac_idx_sample), axis=None)
#test_x = test_x[idx_sample,:]
#y_labels = y_labels[idx_sample]


model = load_model(sys.argv[1])

ps = model.predict(test_x)

right_path = []
wrong_path = []
human = []


print("Either score > 0.7:")
qq = 0.7
te = np.where(ps[:,2] > qq, 2, 1)
seed_labels = np.where(ps[:,0] > qq, 0, te)
cm = confusion_matrix(y_labels, seed_labels)
print(cm)
right_path.append(cm[0,0] + cm[2,2])
wrong_path.append(cm[0,2] + cm[2,0])
human.append(cm[1,0] + cm[1,2])
print(cm[0,0], cm[2,2])
print(cm[0,2], cm[2,0])
print(cm[1,0], cm[1,2])

print("Either score > 0.6:")
qq = 0.6
te = np.where(ps[:,2] > qq, 2, 1)
seed_labels = np.where(ps[:,0] > qq, 0, te)
cm = confusion_matrix(y_labels, seed_labels)
print(cm)
right_path.append(cm[0,0] + cm[2,2])
wrong_path.append(cm[0,2] + cm[2,0])
human.append(cm[1,0] + cm[1,2])

print("Either score > 0.5:")
qq = 0.5
te = np.where(ps[:,2] > qq, 2, 1)
seed_labels = np.where(ps[:,0] > qq, 0, te)
cm = confusion_matrix(y_labels, seed_labels)
print(cm)
right_path.append(cm[0,0] + cm[2,2])
wrong_path.append(cm[0,2] + cm[2,0])
human.append(cm[1,0] + cm[1,2])

print("Either score > 0.46 and largest:")
qq = 0.46
maxes = ps.argmax(axis=1)
#te = np.where(np.logical_and(np.logical_and(ps[:,2] > qq, ps[:,2] > ps[:,1]), ps[:,2] > ps[:,0]), 2, 1)
#seed_labels = np.where(np.logical_and(np.logical_and(ps[:,1] > qq, ps[:,1] > ps[:,2]), ps[:,1] > ps[:,0]), 0, te)
te = np.where(ps[:,2] > qq, 2, 1) #2 if bact > 0.46, else 1
ye = np.where(maxes == 2, te, 1) #2 if bact > 0.46 and max, else 1
ue = np.where(ps[:,0] > qq, 0, 1) #0 if virus > 0.46, else 1
seed_labels = np.where(maxes == 0, ue, ye) #0 if virus > 0.46 & max, 2 if bact > 0.46 and max, else 1
cm = confusion_matrix(y_labels, seed_labels)
print(cm)
right_path.append(cm[0,0] + cm[2,2])
wrong_path.append(cm[0,2] + cm[2,0])
human.append(cm[1,0] + cm[1,2])

print("Human Score < 1/3:")
qq = 1/3
te = np.where(ps[:,0] > ps[:,2], 0, 2)
seed_labels = np.where(ps[:,1] < qq, te, 1)
cm = confusion_matrix(y_labels, seed_labels)
print(cm)
right_path.append(cm[0,0] + cm[2,2])
wrong_path.append(cm[0,2] + cm[2,0])
human.append(cm[1,0] + cm[1,2])

path_sum = [right_path[i] + wrong_path[i] for i in range(len(right_path))]
labels = ["Either > 0.7", "Either > 0.6", "Either > 0.5", "Either > 0.46", "Human < 1/3"]

fig, ax = plt.subplots()
ax.bar(labels, right_path, label='Corret Pathogen', color="g")
ax.bar(labels, wrong_path, label='Opposite Pathogen', color="y", bottom=right_path)
ax.bar(labels, human, label='Human', color="r", bottom=path_sum)

ax.set_ylabel('Total Seeds')
ax.legend()

plt.show()

