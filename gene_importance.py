import pandas as pd
import numpy as np

"""Calculate gene importance for each drug"""

equal_weight = True

def eval_literal(literal):
    """Return ('x', 1) for 'x' and ('x', -1) for '~x'"""
    if literal.startswith('~'):
        return (literal.replace('~','',1), -1)
    else:
        return (literal, 1)

def calc_weights_formula(formula, weights={}):
    """Split formula into literals and calculate weights
       Approach: - each conjunction has weight 1/#conjunctions in formula
                 - each literals has same weight in conjunction
                 - negative propositions (~) have weight -1 * weight of literal
                 - positive propositions have +1 weight * weight of literal
    """
    conjunctions = formula.split('|')
    if equal_weight:
        conjunctions_weight = [(conj.strip(), 1) for conj in conjunctions]
    else:
        conjunctions_weight = [
            (conj.strip(), 1./len(conjunctions)) for conj in conjunctions]

    for conj in conjunctions_weight:
        literals = conj[0].split('&')
        for lit in literals:
            lw = eval_literal(lit.strip())
            if lw[0] not in weights:  # gene has no weight yet
                weights[lw[0]] = 0
            weights[lw[0]] += lw[1] * conj[1]

    return weights


extraLabel = ''
extraLabel = '{}_t0.1'.format(extraLabel)

log_units = True
normalized = False
if normalized:
    extraLabel = '{}_norm'.format(extraLabel)
if log_units:
    extraLabel = '{}_log10'.format(extraLabel)

extraLabel = '{}_limit4'.format(extraLabel)

d = pd.read_csv('ROC_formulae_test{}.csv'.format(extraLabel))
drugs = np.unique(d['Drug'])

saveTable = True

if saveTable:
    if equal_weight:
        extraLabel = '{}_ew'.format(extraLabel)
    outfile = open('gene_importance{}.csv'.format(extraLabel), 'w')
    outfile.write('Drug,num_formula,gene,weight\n')

for drug in drugs:

    print('---- {} ----'.format(drug))
    sel = d['Drug'] == drug
    weights = {}
    for formula in d[sel][' formula']:
        print(formula)
        weights = calc_weights_formula(formula, weights)
        print(weights)

    if saveTable:
        for gene in weights:
            outfile.write(
                '{},{},{},{}\n'.format(
                    drug, np.sum(sel), gene, weights[gene]))

if saveTable:
    outfile.close()
