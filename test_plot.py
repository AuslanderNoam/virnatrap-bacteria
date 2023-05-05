
import sys
import numpy as np
from tensorflow.keras.models import load_model
from sklearn.metrics import confusion_matrix, precision_recall_curve
from tensorflow.keras.metrics import AUC
sys.path.append("/wistar/auslander/Daniel/scripts")
from cnn3class import load_data
import matplotlib.pyplot as plt


test_x, test_y = load_data("test_all_76bp_filt.npz")
int_labels = test_y.argmax(axis=1)
vir_idx = np.where(int_labels == 0)[0]
hum_idx = np.where(int_labels == 1)[0]
bac_idx = np.where(int_labels == 2)[0]
print(vir_idx.shape)
rng = np.random.default_rng(0)
k=1000
vir_idx_sample = rng.choice(vir_idx, size=k, replace=False)
hum_idx_sample = rng.choice(hum_idx, size=k, replace=False)
bac_idx_sample = rng.choice(bac_idx, size=k, replace=False)
idx_sample = np.concatenate((vir_idx_sample, hum_idx_sample, bac_idx_sample), axis=None)
test_x_sample = test_x[idx_sample,:]

model = load_model("final_model_2")

ps = model.predict(test_x_sample)
print("Predicted")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
markers=["^","o","s"]
labels=["Virus", "Human", "Bacteria"]
for i in range(3):
    ax.scatter(ps[i*k:(i+1)*k,0], ps[i*k:(i+1)*k,1], ps[i*k:(i+1)*k,2], marker=markers[i], label=labels[i])

ax.set_xlabel('Virus Score')
ax.set_ylabel('Human Score')
ax.set_zlabel('Bact Score')

plt.show()