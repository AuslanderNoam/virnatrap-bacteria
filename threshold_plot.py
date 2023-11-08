import matplotlib.pyplot as plt

#Based on species-level filter
#sig_GTEx = [90, 70, 62, 59, 54, 53, 50, 45, 41, 39]
#sig_ESCA = [32, 34, 32, 29, 26, 26, 27, 27, 26, 27]
#not_sig =  [23, 14, 12, 10, 10, 11, 9, 9, 11, 10]

#Based on genus-level filter
sig_GTEx = [90,70,65,59,55,54,51,45,44,42]
sig_ESCA = [32,32,31,30,27,26,25,24,24,24]
not_sig =  [23,17,13,10,10,10,11,12,12,12]

sig_sum = [sig_GTEx[i] + sig_ESCA[i] for i in range(len(sig_GTEx))]
labels = [str(i) for i in range(1,11)]

fig, ax = plt.subplots(layout='constrained')
fig.set_size_inches(6.5,4)
rects = ax.bar(labels, sig_GTEx, label='Sig. more prevalent in healthy', color="g")
rects = ax.bar(labels, sig_ESCA, label='Sig. more prevalent in cancer', color="r", bottom=sig_GTEx)
rects = ax.bar(labels, not_sig, label='Not Significant', color="grey", bottom=sig_sum)

ax.set_xlabel('Min. # of Contigs for Detection')
ax.set_ylabel("Number of Genera with >10% Prevalence")
ax.legend(loc='upper right')

#https://www.pythoncharts.com/matplotlib/stacked-bar-charts-labels/
y_offset=-6
for bar in ax.patches:
  ax.text(
      # Put the text in the middle of each bar. get_x returns the start
      # so we add half the width to get to the middle.
      bar.get_x() + bar.get_width() / 2,
      # Vertically, add the height of the bar to the start of the bar,
      # along with the offset.
      bar.get_height() + bar.get_y() + y_offset,
      # This is actual value we'll show.
      round(bar.get_height()),
      # Center the labels and style them a bit.
      ha='center',
      #color='w',
      #weight='bold',
      size=8
  )


plt.savefig("contig_thresh.png", dpi=200)

