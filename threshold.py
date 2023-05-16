import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import matplotlib.pyplot as plt
import sklearn.mixture as skm
import scipy.interpolate
import scipy.optimize

"""Calculate binarization threshold following Knijnenburg et al. 2016"""

plt.ion()

np.random.seed(1)

# all cell lines for which we have genetic data => for LOBICO analysis
drugFilename = 'Supplementary Table 5 LOBICO data for paper.xlsx'

# all cell lines for which we have AUC values & limits => for threshold
AUC_filename = 'AUC_upper_lower_bounds.xlsx'

saveFigure = True
saveTable = True
fontsize = 14
N_upsample = 1000
t = 0.1
log_units = True
normalized = False

extraLabel = '_t{}'.format(t)

if normalized:
    extraLabel = '{}_norm'.format(extraLabel)
if log_units:
    extraLabel = '{}_log10'.format(extraLabel)

AUC = pd.read_excel(AUC_filename, sheet_name='AUC', index_col=0)
AUC_upper = pd.read_excel(AUC_filename, sheet_name='AUC_upper', index_col=0)
AUC_lower = pd.read_excel(AUC_filename, sheet_name='AUC_lower', index_col=0)

drugs = list(AUC.index)
print(f'drugs = {drugs}')

cell_lines = list(AUC)
print(f'cell lines = {cell_lines}')

panel_structure = (3, 2)
n_panels = panel_structure[0] * panel_structure[1]

if saveTable:
    outfile = open('AUC{}.csv'.format(extraLabel), 'w')
    outfile.write('#_drug, num_cell_lines, b, mu, Sigma, theta\n')

    outfile2 = open('Samples{}.csv'.format(extraLabel), 'w')
    outfile2.write('#_drug, num_features, num_cell_lines_thresh, '
                   'num_sensitive_thresh, num_resistant_thresh, '
                   'num_cell_lines_LOBICO, num_sensitive_LOBICO, '
                   'num_resistant_LOBICO\n')

ifig = 0
for idrug, drug in enumerate(drugs):

    dd = pd.read_excel(drugFilename, sheet_name=drug)

    transN = lambda x: x
    if normalized:
        transN = lambda x: x / dd[drug].max()
    transL = lambda x: x
    if log_units:
        transL = lambda x: np.log10(x)
    trans = lambda x: transL(transN(x))

    Samples = dd['Cell']
    AUCval = trans(dd[drug])

    numFeatures = dd.shape[1] - 3
    MM = dd.iloc[:, 3:-1].to_numpy()
    missing = np.zeros(len(Samples), dtype=bool)
    for iSample in range(len(Samples)):
        if (np.any((MM[iSample] != 0) & (MM[iSample] != 1)) or
            np.isnan(AUCval[iSample])):
            missing[iSample] = True
    num_avail_cell_lines = len(Samples) - np.sum(missing)

    if idrug % n_panels == 0:
        fig, axes = plt.subplots(
            panel_structure[0], panel_structure[1], figsize=(9, 9),
            squeeze=True)
        plt.subplots_adjust(wspace=0.3, hspace=0.25)
        ifig += 1
        for axis in axes.flatten():
            axis.axis('off')

    plt.subplot(panel_structure[0], panel_structure[1], 1 + idrug % n_panels)
    plt.gca().axis('on')

    num_cell_lines = np.sum(~np.isnan(AUC.loc[drug, :]))

    print('drug = ', drug, ' has', num_cell_lines, 'valid cell lines')

    # upsampled AUCs for given drug
    up_data = np.zeros((num_cell_lines, N_upsample))

    # Upsampling: we use all cell-lines with AUC values
    icell_line = 0
    for cell_line in cell_lines:

        if np.isnan(AUC.loc[drug, cell_line]):
            continue

        # check that AUC values a consistent for cell lines that are both
        # in drugFile and in AUC_file
        if len(np.where(Samples == cell_line)[0]):
            assert np.isclose(trans(AUC.loc[drug, cell_line]),
                              AUCval[np.where(Samples == cell_line)[0][0]])

        mean = trans(AUC.loc[drug, cell_line])
        lower = trans(AUC_lower.loc[drug, cell_line])
        upper = trans(AUC_upper.loc[drug, cell_line])

        # same as [(upper - mean) + (mean - lower)] * 0.5
        diff = (upper - lower)*0.5
        # should be two sigma given upper & lower are 95% conf. interval limits
        sigma = diff / 2.

        # upsampled AUCs for given drug and cell line
        up_data[icell_line, :] = stats.norm.rvs(mean, sigma, N_upsample)
        icell_line += 1

    # plot histogram of up-sampled data
    plt.hist(up_data.flatten(), bins=80, density=True, color=[0.7, 0.7, 0.7],
             edgecolor=[0.7, 0.7, 0.7])

    # kernel density estimate
    data_min = np.min(up_data)
    data_max = np.max(up_data)
    kernel = stats.gaussian_kde(up_data.flatten(),
                                bw_method=0.5)

    all_sigma = up_data.std()
    xx = np.linspace(data_min-all_sigma, data_max+all_sigma, 1000)
    f = kernel.evaluate(xx)
    f_interp = scipy.interpolate.interp1d(xx, f)
    xxc = 0.5*(xx[1:] + xx[:-1])
    df_interp = scipy.interpolate.interp1d(
        xxc, np.diff(f)/np.diff(xx))
    plt.plot(xx, f, '-b', lw=2)

    mu = xx[np.argmax(f)]
    plt.plot([mu, mu], [0., f_interp(mu)], '-k', lw=3)

    sel = (xx < mu)

    xxc = 0.5*(xx[sel][1:] + xx[sel][:-1])
    xxcc = 0.5*(xxc[1:] + xxc[:-1])
    xxccc = 0.5*(xxcc[1:] + xxcc[:-1])
    df = np.diff(f[sel])/np.diff(xx[sel])
    ddf = np.diff(df) / np.diff(xxc)
    dddf = np.diff(ddf) / np.diff(xxcc)

    # Calculate theta
    found_theta = False

    # 1.
    if not found_theta:
        auxs = np.sign(df[:-1] * df[1:])
        solutions = []
        for iaux, aux in enumerate(auxs):
            theta = xxcc[iaux]
            criterion1a = (np.trapz(f[xx <= theta], xx[xx <= theta]) > 0.05)
            criterion1b = (f_interp(theta) < 0.8*f_interp(mu))
            # print theta,np.trapz(f[xx <= theta],xx[xx <= theta]),criterion1a
            if aux < 0 and criterion1a and criterion1b:
                solutions.append([theta])
        if len(solutions) > 0:
            theta = np.max(solutions)
            print('Criterion 1: found theta =', theta)
            found_theta = True

    # 2.
    if not found_theta:
        auxs = np.sign(ddf[:-1] * ddf[1:])
        solutions = []
        for iaux, aux in enumerate(auxs):
            theta = xxccc[iaux]
            criterion2a = (np.trapz(f[xx <= theta], xx[xx <= theta]) > 0.05)
            criterion2b = (f_interp(theta) < 0.8*f_interp(mu))
            criterion2c = dddf[iaux] > 0
            if aux < 0 and criterion2a and criterion2b and criterion2c > 0:
                solutions.append([theta])
        if len(solutions) > 0:
            theta = np.max(solutions)
            print('Criterion 2: found theta =', theta)
            found_theta = True

    # 3.
    if not found_theta:
        print('Set to theta to', data_min)
        theta = data_min
        found_theta = True

    plt.plot([theta, theta], [0., f_interp(theta)], '-m', lw=3)

    # compute Sigma
    Sigma = np.median(
        np.abs(up_data[(up_data >= theta) & (up_data <= mu)] - mu))

    aux_norm = scipy.stats.norm.pdf(mu, mu, Sigma) / f_interp(mu)
    plt.plot(xx, scipy.stats.norm.pdf(xx, mu, Sigma)/aux_norm, '-k', lw=1.5)
    plt.plot([mu - Sigma, mu],
             1./np.sqrt(2*np.pi*Sigma**2) *
             np.exp(-0.5)/aux_norm*np.array([1, 1]),
             '-k', lw=1.5)

    # Binarization threshold
    cdf = scipy.stats.norm.cdf(xx, mu, Sigma)
    b = np.interp(t, cdf, xx)
    print('Binarization threshold is b =', b)
    plt.plot([b, b], [0., scipy.stats.norm.pdf(b, mu, Sigma)/aux_norm],
             '-r', lw=3)

    num_pos_cell_lines = np.sum(np.logical_and(
        ~np.isnan(AUC.loc[drug, :]), trans(AUC.loc[drug, :]) < b))
    num_neg_cell_lines = np.sum(np.logical_and(
        ~np.isnan(AUC.loc[drug, :]), trans(AUC.loc[drug, :]) >= b))
    print('Total: {} positive {} negative {}'.format(
        num_cell_lines, num_pos_cell_lines, num_neg_cell_lines))

    numSensitive = np.sum(AUCval < b)
    numResistant = np.sum(AUCval >=b)
    if saveTable:
        outfile.write('"{}", {}, {:.5f}, {:.3f}, {:.3f}, {:.3f}\n'.format(
                      drug, num_cell_lines, b, mu, Sigma, theta))
        outfile2.write('{}, {:d}, {:d}, {:d}, {:d}, {:d}, {:d}, {:d}\n'.format(
                      drug, numFeatures,
                      num_cell_lines, num_pos_cell_lines, num_neg_cell_lines,
                      num_avail_cell_lines, numSensitive, numResistant))

    xlabel = 'AUC'
    if normalized:
        xlabel = f'{xlabel} [norm.]'
    if log_units:
        xlabel = rf'$\log_{{10}}$ {xlabel}'
    plt.xlabel(xlabel, fontsize=fontsize-1)

    plt.ylabel('probability density', fontsize=fontsize-1)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)

    plt.tight_layout()

    plt.text(0.03, 0.97, '{} ({} cell lines)'.format(drug, num_cell_lines),
             fontsize=fontsize, verticalalignment='top',
             horizontalalignment='left', transform=plt.gca().transAxes)
    plt.text(0.03, 0.82, 'b={:.4f}'.format(b),
             fontsize=fontsize, verticalalignment='top',
             horizontalalignment='left', transform=plt.gca().transAxes)

    if idrug % n_panels == n_panels-1 or idrug == len(drugs)-1:
        if saveFigure:
            fig.savefig('AUC_f{}{}.pdf'.format(ifig, extraLabel))

plt.show()
if saveTable:
    outfile.close()
    outfile2.close()
