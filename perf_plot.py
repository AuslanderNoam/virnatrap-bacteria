
import sys
import numpy as np
from tensorflow.keras.models import load_model
from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay
#sys.path.append("/wistar/auslander/Daniel/scripts")
from cnn3class import load_data
import matplotlib.pyplot as plt

test_x, test_y = load_data(sys.argv[2])
y_labels = test_y.argmax(axis=1)


model = load_model(sys.argv[1])

ps = model.predict(test_x)

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(6,8), sharex=False)

classes = ["Virus", "Human", "Bacteria"]
for i in range(3):
    PrecisionRecallDisplay.from_predictions(np.where(y_labels == i, 1, 0), ps[:,i], name=f"{classes[i]} vs All", ax=ax1)
ax1.set_xlabel('Recall')
ax1.set_ylabel('Precision')
ax1.set_title("PR Curves")
ax1.set_ylim(0,1)

for i in range(3):
    RocCurveDisplay.from_predictions(np.where(y_labels == i, 1, 0), ps[:,i], name=f"{classes[i]} vs All", ax=ax2)
ax2.set_xlabel('True Positive Rate')
ax2.set_ylabel('False Positive Rate')
ax2.set_title("ROC Curves")
ax2.set_ylim(0,1)

plt.rcParams["font.family"] = "Arial"

plt.rcParams['pdf.fonttype'] = 42

plt.rcParams['axes.facecolor'] = 'none'


fig.tight_layout()
plt.savefig(sys.arv[3])
plt.close()

