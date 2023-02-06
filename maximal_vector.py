"""Basic algorithm solving the maximal vector problem"""
import numpy as np

# See "BEST" algorithm in Godfrey, Shipley & Gryz 2006

def get_maximal_vectors(vectors, gt):
    """find maxima in list l according to 'gt' strict partial order """
    max_vectors = []
    while vectors:
        max_vector = vectors.pop(0)

        vectors_keep = []
        while vectors:  # find a maximum
            vector = vectors.pop(0)
            if not gt(max_vector, vector) and not gt(vector, max_vector):
                vectors_keep.append(vector)
            if gt(vector, max_vector):
                max_vector = vector
        vectors = vectors_keep

        max_vectors.append(max_vector)

        vectors_keep = []
        while vectors:  # clean up
            vector = vectors.pop(0)
            if not gt(max_vector, vector):
                vectors_keep.append(vector)
        vectors = vectors_keep

    return max_vectors

if __name__ == '__main__':

    vectors = [(1, 2), (2, 1), (1, 2), (0, 0.2), (-0.5, 2.5), (1.2, 1.3)]
    def gt(x, y):
        return np.all(np.array(x) >= np.array(y)) and np.any(np.array(x) > np.array(y))
    max_vectors = get_maximal_vectors(vectors, gt)
    print(max_vectors)
    print(set(max_vectors))

    def is_equal(x, y):
        """x and y are 3 vectors corresponding to sens, spec, length of formula"""
        return ((x[0] == y[0]) and (x[1] == y[1]) and (x[2] == y[2]))

    def is_greater(x, y):
        """x and y are 3 vectors corresponding to sens, spec, length of formula"""
        return (((x[0] > y[0]) and (x[1] >= y[1])) or
                ((x[0] >= y[0]) and (x[1] > y[1])) or
                ((x[0] >= y[0]) and (x[1] >= y[1]) and (x[2] < y[2])))

    def is_greater_or_equal(x, y):
        """x and y are 3 vectors corresponding to sens, spec, length of formula"""
        return is_equal(x, y) or is_greater(x, y)

    def get_best_formulae(sens_list, spec_list, form_list):
        """Return indices of entries with optimal sens and spec"""
        n = len(sens_list)
        vector_data =  np.vstack([sens_list, spec_list, form_list, np.arange(n)]).T.tolist()
        print(vector_data)
        max_vectors = get_maximal_vectors(vector_data, is_greater)
        max_vector_indices = np.array(max_vectors)[:, -1].tolist()
        max_vectors = np.array(max_vectors)[:, :-1].tolist()
        return max_vectors, max_vector_indices

    sens_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    spec_list = [0.9, 0.2, 0.7, 0.3, 0.5, 0.5]
    form_list = [1,   2,   3,   2,   1,   2]
    print(get_best_formulae(sens_list, spec_list, form_list))
