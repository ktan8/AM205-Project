# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import networkx as nx
import pandas as pd
import os
import scipy.stats

def sens_prop(adj):
    '''Generate senstitivy of propogration networks as described in 
    Equation 3'''
    adj = np.matrix(adj)
    D1 = np.array(np.sqrt(np.sum(np.fabs(adj), axis = 1)))
    D2 = np.array(np.sqrt(np.sum(np.fabs(adj), axis = 0)))
    
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
    (sd, ud, uu, _) = generate_topology(data)    
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
    dist_ud = dist_pop(create_distance(ud))
    dist_uu = dist_pop(create_distance(uu))
    return(dist_ud, dist_uu)
    
    
def create_distance(mat):
    '''Take adajcent matric and copute all pairs shortest path'''
    np.fill_diagonal(mat, 0)
    g = nx.Graph(mat)
    ds = nx.floyd_warshall_numpy(g)
    return(ds)

def dist_pop(dist_mat):
    '''Distance matrix -> S matrix'''
    return(1 / (1 + dist_mat))
    
def example():
    '''Example generation of data'''
    (data, indexname, ind) = get_data(os.path.join("BIOMD0000000313", "jacobian.csv"))
   # (data, indexname) = get_data(os.path.join("BIOMD0000000404", "jacobian.csv"))

    prop = generate_sense_prop(data)
    dist = generate_distance_influence(data)
    neighbors = generate_topology(data)
    return(prop, dist, (neighbors[0], neighbors[1], neighbors[3]))
    
def compare_mat(mat1, mat2):
    a = np.linalg.norm(mat1 - mat2)
    scipy.stats.spearmanr(mat1, mat2)

