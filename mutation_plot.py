import matplotlib.pyplot as plt
import numpy as np

labels = ["0", "1", "2", "3", "4"]

auprc_dict = {"Bacteria": [0.9069383, 0.89682066, 0.8859337, 0.8741814, 0.8617556],
              "Human": [0.8975821, 0.88160974, 0.8650929, 0.8477688, 0.82959604],
              "Virus": [0.80145264, 0.7792551, 0.7566913, 0.73372316, 0.7109531]}


#Based on https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
fig, ax = plt.subplots(layout='constrained')
fig.set_size_inches(6.5, 4)
x = np.arange(len(labels))  # the label locations
width = 0.3  # the width of the bars
multiplier = 0

for attribute, measurement in auprc_dict.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3, fontsize=8, fmt="%.3f")
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Validation-Set Class One-Versus-All AUPRC')
ax.set_xlabel("Number of Simulated Mutations")
ax.set_xticks(x + width, labels)
ax.legend(loc='upper left', ncol=3)
ax.set_ylim(0.5, 1.0)

plt.savefig("mut_auprc.png", dpi=200)
