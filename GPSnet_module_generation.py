# -*- coding: utf-8 -*-
'''
Description: random module generate by GPSnet algorithm.
Input 1: node_score.txt:  column 1: protein entrez id, column 2: node score -- in our case, |log2FC|
Input 2: pickle file: ppi.p which stores the protein-protein interactome, here we store it in a dictionary, with the protein as the key and their first-order neighbor as the value
'''

import os
import time
import random
import pickle
import numpy as np
from math import sqrt
from scipy.stats import hypergeom

P_CUTOFF = 1E-2


DIR_DATA = '../data'
DIR_MODULE = '../modules'


with open(os.path.join(DIR_DATA, 'node_score.txt')) as fi:
    next(fi)
    NODE_ID = []
    NODE_SCORE = []
    for line in fi:
        id, score = line.strip("\n").split("\t")
        NODE_ID.append(id)
        NODE_SCORE.append(score)

NODE_ID = np.array(NODE_ID, dtype=np.int)  # Array: pos -> id
NODE_SCORE = np.array(NODE_SCORE, dtype=np.float)  # Array: pos -> score

with open(os.path.join(DIR_DATA, 'ppi.p'), 'rb') as fi:
    NODE_NEIGHBOR = pickle.load(fi)  # ['p1':{n1, n2, ...}, ...]


NODE_DEGREE = np.array([len(NODE_NEIGHBOR[node]) for node in NODE_ID])  # Array: pos -> degree


N_NODES = len(NODE_ID)
MEAN_NODE_SCORE = NODE_SCORE.mean()
ID2POS = {j: i for i, j in enumerate(NODE_ID)}  # {id -> pos}
NODE_ID_SET = set(NODE_ID)




# ----------------------------------------------------------
def P_Connectivity(k_m, k_i, n):
    p_rv = hypergeom(N_NODES, n, k_i)
    return p_rv.cdf(min(k_i, n)) - p_rv.cdf(k_m - 1)

def FindModule(idx):
# def FindModule(idx, seed):
    # random.seed(seed)

    last_added_pos = random.randrange(0, N_NODES)
    module = [last_added_pos]  # [pos1, pos2, ...]
    module_id = {NODE_ID[last_added_pos]}  # {id1, id2, ...}
    score_sum = 0
    neighbors_id = set()


    while True:
        print(".", end="")

        # -------- add note here -------- #
        module_size = len(module)
        score_sum += NODE_SCORE[last_added_pos]
        module_score = (score_sum - module_size * MEAN_NODE_SCORE) / sqrt(module_size)

        neighbors_id |= NODE_NEIGHBOR[NODE_ID[last_added_pos]]
        neighbors_id -= module_id
        if not len(neighbors_id):
            break
        neighbors_pos = np.array([ID2POS[id] for id in sorted(neighbors_id)])  ###

        # -------- select nodes that could increase module score -------- #
        qualified_pos = neighbors_pos[module_score * sqrt(module_size) + NODE_SCORE[neighbors_pos] > module_score * sqrt(module_size + 1) + MEAN_NODE_SCORE]
        if not len(qualified_pos):
            break

        # -------- select nodes whose connectivity significance p-value <= P_CUTOFF -------- #
        k_m_ls = [len(NODE_NEIGHBOR[NODE_ID[pos]] & module_id) for pos in qualified_pos]
        k_i_ls = NODE_DEGREE[qualified_pos]
        qualified_csp_pos = [pos for pos, k_m, k_i in zip(qualified_pos, k_m_ls, k_i_ls) if P_Connectivity(k_m, k_i, module_size) <= P_CUTOFF]
        if not len(qualified_csp_pos):
            break

        # -------- add note here -------- #
        qualified_csp_score = NODE_SCORE[qualified_csp_pos]
        max_pos_index = np.where((qualified_csp_score == np.max(qualified_csp_score)))[0]
        last_added_pos = qualified_csp_pos[random.choice(max_pos_index)]
        module.append(last_added_pos)
        module_id.add(NODE_ID[last_added_pos])

    print("")

    with open(os.path.join(DIR_MODULE, 'module_%s.txt' % idx), "w") as fo:
        fo.write("node_id\n")
        fo.write("\n".join(str(i) for i in NODE_ID[module]))


# ----------------------------------------------------------
for idx in range(0, 100000):
    start_time = time.time()
    # FindModule(idx, seed=idx)
    FindModule(idx)
    print("--- Experiment %s takes %.1f seconds ---" % (idx, time.time() - start_time))
