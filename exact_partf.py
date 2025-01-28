from tree_scorer import calc_prob, pf_cond_on_one_tree
from utils import all_trees, newick_to_subtrees, binary_trees
import pandas as pd
from decimal import Decimal
import numpy as np
import argparse


def exact_partf(P, cells, mutation, df):

    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    cond_m = np.where(df.columns == mutation)[0][0]

    n = len(df.index)

    numerator = Decimal(0)
    denominator = Decimal(0)
    for tree in all_trees(n):
        sts = newick_to_subtrees(tree)
        nontrivial = sts[n + 1 : -1]
        order = len(nontrivial)
        p = calc_prob(P, sts, order)

        cond = pf_cond_on_one_tree(P, sts, cond_c, cond_m)
        numerator +=  Decimal(cond[0] / cond[1]) * Decimal(p)
        denominator += Decimal(p)

    return (numerator / denominator)

def our_partf(P, cells, mutation, df):

    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    cond_m = np.where(df.columns == mutation)[0][0]

    n = len(df.index)

    numerator = Decimal(0)
    denominator = Decimal(0)
    for tree in binary_trees(n):
        sts = newick_to_subtrees(tree)
        order = 0
        p = calc_prob(P, sts, order)

        cond = pf_cond_on_one_tree(P, sts, cond_c, cond_m)
        numerator +=  Decimal(cond[0] / cond[1]) * Decimal(p)
        denominator += Decimal(p)

    return (numerator / denominator)

def main(args):

    
    path = args.input_matrix
    cells = args.cells
    mutation = args.mutation
    alpha = args.alpha
    beta = args.beta

    # path = "/Users/john/Desktop/PartF-Journal-Jan_18_2024/simNo_1-n_8-m_8-fp_0.001-fn_0.1-na_0.SC.after_noise"
    # alpha=0.001
    # beta=0.1
    # mutation = "mut2"
    # cells = ["cell_1", "cell_2", "cell_3"]



    df = pd.read_csv(path, sep="\t", index_col=[0]).sort_values(by=["cell_id_x_mut_id"])
    idx_to_cells = df.index
    cells_to_idx = {idx_to_cells[i]:i for i in range(len(idx_to_cells))}
    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    pstr = None

    if args.exact:
        pstr = str(exact_partf(P, cells, mutation, df))
    else:
        pstr = str(our_partf(P, cells, mutation, df))

    output = [args.path, pstr, args.exact]
    print("\t".join(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run.py')

    parser.add_argument("-i", "--input_matrix", type=str,                                                        
                        help="Path to input genotype matrix where rows correspond to cells/sublines and columns correspond to mutations. See repo examples for formatting.", required=True)
    parser.add_argument("-c", "--cells", type=str,                                                        
                        help="List of cells (comma separated)", required=True)
    parser.add_argument("-m", "--mutation", type=str,                                                        
                        help="Name of the mutation (column) in the matrix", required=True)
    parser.add_argument("-fp", "--alpha", type=float,                                                        
                        help="False-positive rate (alpha in the paper)", required=True)
    parser.add_argument("-fn", "--beta", type=float,                                                        
                        help="False-negative rate (beta in the paper)", required=True)
    parser.add_argument("-e", "--exact", type=bool,                                                        
                        help="True for exact, False for binary", required=True)
    
    
    main(parser.parse_args())
