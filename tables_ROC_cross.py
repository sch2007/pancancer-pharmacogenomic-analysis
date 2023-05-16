import pandas as pd
import numpy as np

from maximal_vector import get_maximal_vectors

np.random.seed(1)

extraLabel = ''

extraLabel = '{}_t0.1'.format(extraLabel)
niter = 5
nfold = 4

log_units = True
normalized = False
if normalized:
    extraLabel = '{}_norm'.format(extraLabel)
if log_units:
    extraLabel = '{}_log10'.format(extraLabel)

KMs = [[1, 1], [1, 2], [2, 1], [1, 3], [3, 1], [2, 2], [1, 4], [4, 1]]
extraLabel = '{}_limit4'.format(extraLabel)

KM_max = np.max(np.array(KMs))
d1 = pd.read_csv('Validation_n{}_k{}_KM{}{}.csv'.format(niter, nfold, KM_max, extraLabel))
d1['sens_train'] =  d1['num_true_pos_train'] / d1['num_pos_train']
d1['spec_train'] =  d1['num_true_neg_train'] / d1['num_neg_train']
d1['sens_val'] =  d1['num_true_pos_val'] / d1['num_pos_val']
d1['spec_val'] =  d1['num_true_neg_val'] / d1['num_neg_val']

d2 = pd.read_csv('Testing_n{}_KM{}{}.csv'.format(niter, KM_max, extraLabel))
d2['sens_trainval'] =  d2['num_true_pos_trainval'] / d2['num_pos_trainval']
d2['spec_trainval'] =  d2['num_true_neg_trainval'] / d2['num_neg_trainval']
d2['sens_test'] =  d2['num_true_pos_test'] / d2['num_pos_test']
d2['spec_test'] =  d2['num_true_neg_test'] / d2['num_neg_test']

drugs = np.unique(d2['drug'])

saveTable = True

# compute cumulative maximum of list
# or n-th previous cumulative maximum
def cummax(l, n=0):
    ret = -np.infty*np.ones(len(l))
    for i in range(len(l)):
        if i+1-n > 0:
            ret[i] = np.max(l[:i+1-n])
    return ret

def get_best_formulae(sens_list, spec_list, form_list):
    """Return indices of entries with optimal sens and spec"""

    def is_greater(x, y):
        return (((x[0] > y[0]) and (x[1] >= y[1])) or
                ((x[0] >= y[0]) and (x[1] > y[1])) or
                ((x[0] >= y[0]) and (x[1] >= y[1]) and (x[2] < y[2])))

    n = len(sens_list)
    vector_data =  np.vstack(
        [sens_list, spec_list, form_list, np.arange(n)]).T.tolist()
    max_vectors = get_maximal_vectors(vector_data, is_greater)
    return max_vectors

def calc_ROC(xx, max_vectors):

    def sort_key(x):
        return x[1]  # sort by specificity

    max_vectors_sorted = sorted(max_vectors, key=sort_key, reverse=True)
    old_sens = 0.
    old_invspec = 0
    yy = 0.*xx
    for max_vector in [(0, 1, 0, 0)] + max_vectors_sorted + [(1, 0, 0, 0)]:
        sens = max_vector[0]
        invspec = 1-max_vector[1]
        if (invspec > old_invspec):
            sel = (xx >= old_invspec) & (xx < invspec)
            yy[np.argmin(np.abs(xx-invspec))] = sens
            yy[sel] = old_sens
        old_sens = sens
        old_invspec = invspec
    yy = np.maximum(yy, xx)
    yy[0] = 0

    return yy

def get_AUC(xx, yy):
    return np.trapz(yy, xx)

if saveTable:
    outfile = open('ROC_AUC{}.csv'.format(extraLabel), 'w')
    outfile.write('Drug, AUC_test, e_AUC_test, AUC_train, e_AUC_train\n')

    outfile2 = open('ROC_formulae_train{}.csv'.format(extraLabel), 'w')
    outfile2.write('Drug, iteration, sens_train, spec_train, iformula, formula\n')

    outfile3 = open('ROC_formulae_test{}.csv'.format(extraLabel), 'w')
    outfile3.write('Drug, iteration, sens_test, spec_test, iformula, formula\n')

for drug in drugs:

    print('---- {} ----'.format(drug))

    N_interp = 1000
    xx = np.linspace(0, 1, N_interp)
    yy_trainval = np.zeros((niter, len(xx)))
    yy_test = np.zeros_like(yy_trainval)

    best_vector_indices_trainval = []
    best_vector_forms_trainval = set()
    KMs_opt = np.zeros((niter, 2), dtype=int)
    for iiter in range(1, niter+1):

        ## ---- Find optimal hyper parameters
        AUC_val_avg = np.zeros(len(KMs))
        AUC_train_avg = np.zeros(len(KMs))
        for iKM, KM in enumerate(KMs):

            AUC_val = np.zeros(nfold)
            AUC_train = np.zeros(nfold)
            for ifold in range(1, nfold+1):

                # 1) load Sensitivity and Specificity of training and testing
                sel = ((d1['drug'] == drug) &
                       (d1['niter'] == niter) & (d1['iiter'] == iiter) &
                       (d1['nfold'] == nfold) & (d1['ifold'] == ifold) &
                       (~np.isnan(d1['sens_val'])) & (~np.isnan(d1['spec_val'])) &
                       (d1['clauses'] == KM[0]) & (d1['literals'] == KM[1]))

                if np.sum(sel) == 0:
                    AUC_val[ifold-1] = np.nan
                    AUC_train[ifold-1] = np.nan
                    continue

                dd = d1[sel][
                    ['clauses', 'literals',
                     'sens_train', 'spec_train',
                     'sens_val', 'spec_val', 'formula']]
                dd = dd.drop_duplicates()

                dd['len'] = [len(x.replace('&', ' ').replace('|', ' ').split()) for x in dd['formula'].to_list()]

                # 2) identify the "best" formulae for the training data
                max_vectors_train = get_best_formulae(
                    dd['sens_train'].to_numpy(), dd['spec_train'].to_numpy(),
                    (dd['len']).to_numpy())
                max_vector_indices = np.array(max_vectors_train)[:, -1].astype(int).tolist()

                max_vectors_val = get_best_formulae(
                    dd.iloc[max_vector_indices]['sens_val'].to_numpy(),
                    dd.iloc[max_vector_indices]['spec_val'].to_numpy(),
                    (dd.iloc[max_vector_indices]['len']).to_numpy())
                max_vector_val_indices = np.array(max_vectors_val)[:, -1].astype(int).tolist()

                # 3) calculate ROC & AUC for validation data
                yy_val  = calc_ROC(xx, max_vectors_val)
                AUC_val[ifold-1] = get_AUC(xx, yy_val)

                yy_train  = calc_ROC(xx, max_vectors_train)
                AUC_train[ifold-1] = get_AUC(xx, yy_train)

            AUC_val_avg[iKM] = np.nanmean(AUC_val)
            AUC_train_avg[iKM] = np.nanmean(AUC_train)

        ## ---- Apply model to test data with optimal hyper parameters

        # 4) pick K&M values with best AUC_val
        print('-- Iteration = {} --'.format(iiter))
        print('training(AUC) = ' +
              ', '.join(list(map('{:.3g}'.format, AUC_train_avg))))
        print('validation(AUC) = ' +
              ', '.join(list(map('{:.3g}'.format, AUC_val_avg))))

        iKM_opt = np.argmax(AUC_val_avg)
        KMs_opt[iiter-1] = KMs[iKM_opt]
        K_opt = KMs[iKM_opt][0]
        M_opt = KMs[iKM_opt][1]
        print('highest val(AUC) K={} M={} w/ val(AUC)={:.3g} train(AUC)={:.3g}'.format(
            K_opt, M_opt, AUC_val_avg[iKM_opt], AUC_train_avg[iKM_opt]))

        # 5) load Sensitivity and Specificity of training and testing
        sel = ((d2['drug'] == drug) &
               (d2['niter'] == niter) & (d2['iiter'] == iiter) &
               (~np.isnan(d2['sens_test'])) & (~np.isnan(d2['spec_test'])) &
               (d2['clauses'] == K_opt) & (d2['literals'] == M_opt))

        if np.sum(sel) == 0:
            yy_trainval[iiter-1, :] = np.nan
            yy_test[iiter-1, :] = np.nan
            continue

        dd = d2[sel][
            ['clauses', 'literals',
             'sens_trainval', 'spec_trainval',
             'sens_test', 'spec_test', 'formula']]
        dd = dd.drop_duplicates()

        # 6) identify the "best" formulae for the training+validation data
        dd['len'] = [len(x.replace('&', ' ').replace('|', ' ').split()) for x in dd['formula'].to_list()]

        max_vectors_trainval = get_best_formulae(
            dd['sens_trainval'].to_numpy(), dd['spec_trainval'].to_numpy(),
            (dd['len']).to_numpy())
        max_vector_indices = np.array(max_vectors_trainval)[:, -1].astype(int).tolist()

        # store line index of all "best" formula based on train+val
        best_vector_indices_trainval.extend(dd.index[max_vector_indices].to_numpy())
        best_vector_forms_trainval.update(dd.iloc[max_vector_indices]['formula'].tolist())

        max_vectors_test = get_best_formulae(
            dd.iloc[max_vector_indices]['sens_test'].to_numpy(),
            dd.iloc[max_vector_indices]['spec_test'].to_numpy(),
            dd.iloc[max_vector_indices]['len'].to_numpy())
        max_vector_test_indices = np.array(max_vector_indices)[
            np.array(max_vectors_test)[:, -1].astype(int).tolist()]

        # 7) calculate ROC for train & test data
        yy_trainval[iiter-1, :]  = calc_ROC(xx, max_vectors_trainval)
        yy_test[iiter-1, :]  = calc_ROC(xx, max_vectors_test)

        print('testing(AUC)={:.3g} '.format(
            get_AUC(xx, yy_test[iiter-1])), end='')

        print('train+val(AUC)={:.3g}'.format(
            get_AUC(xx, yy_trainval[iiter-1])))

        if saveTable:
            for ii, i in enumerate(max_vector_indices):
                outfile2.write('{}, {:d}, {:.3f}, {:.3f}, {:d}, {}\n'.format(
                    drug, iiter, dd.iloc[i]['sens_trainval'],
                    dd.iloc[i]['spec_trainval'], ii, dd.iloc[i]['formula']))
            for ii, i in enumerate(max_vector_test_indices):
                outfile3.write('{}, {:d}, {:.3f}, {:.3f}, {:d}, {}\n'.format(
                    drug, iiter, dd.iloc[i]['sens_test'],
                    dd.iloc[i]['spec_test'], ii, dd.iloc[i]['formula']))

    # 8) average ROC
    yy_trainval_avg = np.nanmean(yy_trainval, axis=0)
    yy_test_avg = np.nanmean(yy_test, axis=0)

    AUC_trainval = np.zeros(niter)
    AUC_test = np.zeros(niter)
    for iiter in range(niter):
        AUC_trainval[iiter] = get_AUC(xx, yy_trainval[iiter])
        AUC_test[iiter] = get_AUC(xx, yy_test[iiter])

    print('Train+Val: AUC = {:.3g} +- {:.3g}'.format(
        np.nanmean(AUC_trainval),np.nanstd(AUC_trainval)))
    print('Test: AUC = {:.3g} +- {:.3g}'.format(
        np.nanmean(AUC_test),np.nanstd(AUC_test)))

    if saveTable:
        outfile.write('{}, {:.3f}, {:.3f}, {:.3f}, {:.3f}\n'.format(
                      drug, np.nanmean(AUC_test), np.nanstd(AUC_test),
                      np.nanmean(AUC_trainval), np.nanstd(AUC_trainval)))
        outfile.flush()

if saveTable:
    outfile.close()
    outfile2.close()
    outfile3.close()
