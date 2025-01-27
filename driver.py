import pandas as pd
from exact_partf import our_partf

path = "/Users/john/Desktop/PartF-Journal-Jan_18_2024/simNo_1-n_8-m_8-fp_0.001-fn_0.1-na_0.SC.after_noise"
alpha=0.001
beta=0.1
mutation = "mut2"
cells = ["cell_1", "cell_2", "cell_3"]



df = pd.read_csv(path, sep="\t", index_col=[0]).sort_values(by=["cell_id_x_mut_id"])
idx_to_cells = df.index
cells_to_idx = {idx_to_cells[i]:i for i in range(len(idx_to_cells))}
I_mtr = df.values
t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
P = t1 + t2
P[I_mtr == 3] = 0.5

print("Ours should converge to: " + str(our_partf(P, cells, mutation, df)))