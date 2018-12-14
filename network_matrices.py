# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import networkx as nx
import pandas as pd
import os
import scipy.stats
import matplotlib.pyplot as plt

def sens_prop(adj):
    '''Generate senstitivy of propogration networks as described in 
    Equation 3'''
    adj = np.matrix(adj)
    D1 = np.array(np.sqrt(np.sum(np.abs(adj), axis = 1)))
    D2 = np.array(np.sqrt(np.sum(np.abs(adj), axis = 0)))
    
    def root_s(x):
        return(np.where(x == 0, 0, 1/x))
    D1 = root_s(D1)
    D2 = root_s(D2)
    D1 = np.diagflat(D1)
    D2 = np.diagflat(D2)
    Wp = D1* adj * D2
    alpha = 0.9
    S = (1 - alpha) * np.linalg.inv((np.eye(len(D1)) - alpha * Wp))
    return(S)

def senseBio(J):
    I = np.eye(J.shape[0])
    IJ = np.linalg.inv(I - J)
    for i in range(0, len(IJ)):
        IJ[:, i] = IJ[:, i] / IJ[i, i]
    return(IJ)

def get_data(path):
    '''Read jacobian/model data'''
    data =  pd.read_csv(path, index_col = 0)
    mat = data.values
    to_check = np.sum(np.fabs(mat), axis = 0)
    to_check = to_check > 0
    ind = np.where(to_check)
    #print(ind)
    indicies = data.columns
    data = data[indicies[ind]].iloc[ind]

    res_mat = data.values
    res_indicies = data.columns
    #mat = data.values
    #mat = mat.T
    res_mat = res_mat.T
    return(res_mat, res_indicies, data)

def generate_topology(data):
    '''For a given adjaceny matrix - convert to a simple -1, 0, 1 topology'''
    signed_directed = np.sign(data)
    unsigned_directed = np.abs(np.sign(data))
    unsigned_undirected = np.sign(unsigned_directed + unsigned_directed.T)
    unsigned_undirected_t = (unsigned_directed + unsigned_directed.T)

    return(signed_directed, unsigned_directed, unsigned_undirected, unsigned_undirected_t)

def generate_sense_prop(data):
    '''Take topolgoy and generate_propagation S'''
    (sd, ud, uu, __) = generate_topology(data)    
    sens_sd = sens_prop(sd)
    sens_ud = sens_prop(ud)
    sens_uu = sens_prop(uu)
    return(sens_sd, sens_ud, sens_uu)

def convertdf(mat, indicies):
    s3 = pd.DataFrame(mat).set_index(indicies)
    s3.columns = indicies
    return(s3)

def generate_distance_influence(data):
    '''Take topolgy and generate distance_influence'''
    (_, ud, uu, _) = generate_topology(data)
    dist_ud = dist_pop(create_distance(ud.T)).T
    dist_uu = dist_pop(create_distance(uu.T)).T
    return(dist_ud, dist_uu)
    
    
def create_distance(mat):
    '''Take adajcent matric and copute all pairs shortest path'''
    np.fill_diagonal(mat, 0)
    g = nx.DiGraph(mat)
    ds = nx.floyd_warshall_numpy(g)
    return(ds)

def dist_pop(dist_mat):
    '''Distance matrix -> S matrix'''
    return(1 / (1 + dist_mat))
    
def example():
    '''Example generation of data'''
    path = os.path.join("BIOMD0000000313", "jacobian.csv")
    return(run_many(path))
    
def run_many(path):
    '''Run throw the exisitng models'''
    (data, indexname, df) = get_data(path)
    prop = generate_sense_prop(data)
    dist = generate_distance_influence(data)
    neighbors = generate_topology(data)
    return(prop, dist, (neighbors[0], neighbors[1], neighbors[3]))

def compare_mat(mat1, mat2):
    a = np.linalg.norm(mat1 - mat2)
    b = scipy.stats.spearmanr(np.abs(mat1), np.abs(mat2), axis = None)
    c= scipy.stats.spearmanr((mat1), (mat2), axis = None)

    return((a, b[0], c[0]))


def example2():
    (data, indexname, ind) = get_data(os.path.join("BIOMD0000000313", "jacobian.csv"))
    sensB = senseBio(data)
    (prop, dis, neighbors )= example()
    off = lambda : np.random.random() * 1e-7
    (cor_prob_sd, cor_prop_ud, cor_prob_uu) = (compare_mat(sensB, prop[0] ), compare_mat(sensB, prop[1] + off()),  compare_mat(sensB, prop[2] + off()))
    (cor_dist_ud, cor_dist_uu) = (compare_mat(sensB, dis[0] + off()), compare_mat(sensB, dis[1] + off()))
    (cor_neigh_sd, cor_neigh_ud, cor_neigh_uu) = (compare_mat(sensB, neighbors[0].T ), compare_mat(sensB, neighbors[1].T) , compare_mat(sensB, neighbors[2].T ))    
    
    (eign0, eign1, eign2) = (np.linalg.eig(neighbors[0].T)[1], np.linalg.eig(neighbors[1].T)[1], np.linalg.eig(neighbors[2].T)[1])
    D1 = np.sum(neighbors[2].T, axis = 1)
        
    D2 = np.sum(neighbors[1].T, axis = 1)
    L1 = neighbors[2].T - D1
    l1 = np.linalg.eig(L1)[1]
    L2 = neighbors[1].T - D2
    l2 = np.linalg.eig(L2)[1]
    

    (corrE0, corrE1, corrE2) = (compare_mat(sensB, eign0), compare_mat(sensB, eign1), compare_mat(sensB, eign2))
    (lapE1, lapE2) = (compare_mat(sensB, l1), compare_mat(sensB, l2))
    (lap1, lap2) = (compare_mat(sensB, L1), compare_mat(sensB, L2))
    #  D1 = (np.sum(np.abs(), axis = 1))
   
    objects = ('Prop sign', 'Propogation_directed_unsigned',
               'Propogation_undirected', 'Distance_directed', 'Distance_undirected', 
               'Neighbor directed', 'Neighbor undirected', 'total neighbor', 
               'Laplacian sign', 'Laplacian unsigned',
               'Eigenvector sign', 'Eigenvector directed', 'eigenvector undirected', \
               'Laplacian_eigen dir', 'Lapaclain_eign_undir')
    y_pos = np.arange(len(objects))
    diff_data = [cor_prob_sd, cor_prop_ud, cor_prob_uu, cor_dist_ud, cor_dist_uu, cor_neigh_sd, cor_neigh_ud, cor_neigh_uu, lap1, lap2, corrE0, corrE1, corrE2, lapE1, lapE2]
    corrs = [x[1] for x in diff_data]
    plt.bar(y_pos, corrs, align='center', alpha=0.5)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Correlation with biomodel')
    plt.show()

    corrs2 = [x[2] for x in diff_data]
    plt.bar(y_pos, corrs2, align='center', alpha=0.5)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Signed Correlation with biomodel')
    plt.show()
    
    #Calculate the number of equal signed parts. 
    q = lambda inp : np.sum(np.sign(inp) == np.sign(sensB)) / (sensB.shape[0] * sensB.shape[1])
    sign_match = [q(prop[0]), q(prop[1]), q(prop[2]), q(dis[0]), q(dis[1]), 
                  q(neighbors[0]), q(neighbors[1]), q(neighbors[2]), q(L1), 
                  q(L2), q(eign0), q(eign1), q(eign2), q(l1), q(l2)]
    

    plt.bar(y_pos, sign_match, align='center', alpha=0.5)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Percent agrement in direction')
    plt.show()
    

    #Norm distance
    
(_, v) = np.linalg.eig(neighbors[0])
