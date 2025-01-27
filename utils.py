import treeswift as ts
import numpy as np
import pandas as pd
from phylo2vec import vecToNewick, phylovecdomain
from cf_mat_to_newick import cf_to_newick
from itertools import chain, combinations

def is_df_in_list(df, df_list):
    for df_item in df_list:
        if df.equals(df_item):
            return True
    return False

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def consistent_subtrees(nwk, n):
    sts = newick_to_subtrees(nwk)
    trivial = sts[0:n + 1] + [sts[-1]]
    nontrivial = sts[n + 1:-1]
    pset = list(powerset(nontrivial))

    ret = []
    for p in pset:
        mat = trivial[:-1] + list(p) + [trivial[-1]]
        ret += [pd.DataFrame(mat).transpose()]
    return ret

def binary_trees(n):
    vecs = phylovecdomain(n)
    newicks = list(map(lambda x: vecToNewick(x) + ";", vecs))

    return newicks

def all_trees(n):
    binary = binary_trees(n)
    ls = []
    for tree in binary:
        ls += [tree]
        for cf_mat in consistent_subtrees(tree, n)[:-1]:
            # Sort cf_mat to be safe?
            new_tree = cf_to_newick(cf_mat)
            if new_tree not in ls:
                ls += [new_tree]
    return(ls)


def newick_to_subtrees(newick):
    # input: newick string of tree with n leaves labelled 0,1,2,...,n-1

    tree = ts.read_tree_newick(newick)
    leaves = tree.labels(leaves=True, internal=False)
    n = len(list(leaves))

    subtrees = []
    for node in tree.traverse_preorder():
        subtree = np.zeros(n, dtype=np.int8)
        leaves = [int(leaf.get_label()) for leaf in node.traverse_leaves()]
        subtree[leaves] = 1
        subtrees += [subtree]

    subtrees += [np.zeros(n, dtype=np.int8)]
    subtrees = [np.array(x) for x in {(tuple(e)) for e in subtrees}]
    subtrees.sort(key=lambda x : (sum(x), int(''.join(map(str,x)), 2)))
    return subtrees