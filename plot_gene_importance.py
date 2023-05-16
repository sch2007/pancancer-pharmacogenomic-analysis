import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""Plot Figure 5C """

plt.ion()
fontsize=10.5

equal_weight = True
extraLabel = ''
extraLabel = '{}_t0.1'.format(extraLabel)

log_units = True
normalized = False
if normalized:
    extraLabel = '{}_norm'.format(extraLabel)
if log_units:
    extraLabel = '{}_log10'.format(extraLabel)

extraLabel = '{}_limit4'.format(extraLabel)
if equal_weight:
    extraLabel = '{}_ew'.format(extraLabel)

d = pd.read_csv('gene_importance{}.csv'.format(extraLabel))

drugs = np.unique(d['Drug'])

cell = data = pd.read_excel(
    ('Supplementary Table 5 LOBICO data for paper.xlsx'), index_col=0)

genes = cell.keys()[2:]

arr = np.zeros((len(drugs), len(genes)))

for i, drug in enumerate(drugs):
    for j, gene in enumerate(genes):
        sel = (d['Drug'] == drug) & (d['gene'] == gene)
        assert (np.sum(sel) == 0) | (np.sum(sel) == 1)
        if np.sum(sel) == 1:
            arr[i, j] = d['weight'][sel] / d['num_formula'][sel]

fig = plt.figure(figsize=(10, 4))
ax = plt.gca()
im = ax.imshow(arr, cmap='bwr')

 # Create colorbar
cbar = ax.figure.colorbar(im, ax=ax, shrink=0.75)
cbar.ax.set_ylabel('weight', rotation=-90, va="bottom")

ax.set_xticks(np.arange(len(genes)))
ax.set_yticks(np.arange(len(drugs)))
ax.set_xticklabels(genes, fontsize=fontsize-2)
ax.set_yticklabels(drugs)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# We want to show all ticks...
ax.set_xticks(np.arange(len(genes)+1)-0.5, minor=True)
ax.set_yticks(np.arange(len(drugs)+1)-0.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=1.5)
ax.tick_params(which="minor", bottom=False, left=False)

im.set_clim([-1, 1])

#ax.set_title('Relative occurrence of genes in LOBICO formulae')

fig.tight_layout()
plt.savefig('gene_importance{}.pdf'.format(extraLabel))
