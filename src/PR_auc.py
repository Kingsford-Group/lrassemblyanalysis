# PR_auc.py
#
# Laura Tung
#
# Usage: PR_auc.py <roc_filename>

import sys
import numpy as np
from sklearn.metrics import auc


def get_PR_AUC(loaded_roc):

    PR_auc = auc(loaded_roc[:, 0], loaded_roc[:, 1])/10000

    return PR_auc


def load_data(dataset):

    loaded_roc = np.loadtxt(dataset, usecols=(12, 15))

    return loaded_roc


if __name__ == "__main__":
    
    loaded_roc = load_data(sys.argv[1])

    PR_auc = get_PR_AUC(loaded_roc)
    print("PR_AUC = "+str(PR_auc))
