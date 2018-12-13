# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import networkx as nx
import pandas as pd
import os

def sens_prop(adj):
    '''Generate senstitivy of propogration networks as described in 
    Equation 3'''
    adj = np.matrix(adj)
    D1 = np.sqrt(np.sum(np.fabs(adj), axis = 1))
    D2 = np.sqrt(np.sum(np.fabs(adj), axis = 0))
    D1 = np.diagflat(D1)
    D2 = np.diagflat(D2)
    Wp = D1* adj * D2
    alpha = 0.9
    S = (1 - alpha) * np.linalg.inv((np.eye(len(D1)) - alpha * Wp))
    return(S)
    
def get_data(path):
    '''Read jacobian/model data'''
    data =  pd.read_csv(path, index_col = 0)
    indicies = data.columns
    mat = data.values
    return(mat, indicies)

def generate_topology(data):
    '''For a given adjaceny matrix - convert to a simple -1, 0, 1 topology'''
    signed_directed = np.sign(data)
    unsigned_directed = np.abs(np.sign(data))
    unsigned_undirected = np.sign(unsigned_directed + unsigned_directed.T)
    return(signed_directed, unsigned_directed, unsigned_undirected)

def generate_sense_prop(data):
    '''Take topolgoy and generate_propagation S'''
    (sd, ud, uu) = generate_topology(data)    
    sens_sd = sens_prop(sd)
    sens_ud = sens_prop(ud)
    sens_uu = sens_prop(uu)
    return(sens_sd, sens_ud, sens_uu)

def generate_distance_influence(data):
    '''Take topolgy and generate distance_influence'''
    (_, ud, uu) = generate_topology(data)
    dist_ud = dist_pop(create_distance(ud))
    dist_uu = dist_pop(create_distance(uu))
    return(dist_ud, dist_uu)
    
    
def create_distance(mat):
    '''Take adajcent matric and copute all pairs shortest path'''
    mat.fill_diagonal(0)
    g = nx.Graph(mat)
    ds = nx.floyd_warshall_numpy(g)
    return(ds)

def dist_pop(dist_mat):
    '''Distance matrix -> S matrix'''
    return(1 / (1 + dist_mat))
    
def example():
    '''Example generation of data'''
    (data, indexname) = get_data(os.path.join("BIOMD0000000313", "jacobian.csv"))
    prop = generate_sense_prop(data)
    dist = generate_distance_influence(data)
    neighbors = generate_topology(data)
    return(prop, dist, neighbors)
