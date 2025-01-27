from utils import consistent_subtrees, all_trees
from cf_mat_to_newick import cf_to_newick

for i in range(2,7):
    print(len(all_trees(i)))

