# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import networkx as nx
import pandas as pd
import os
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
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
    '''Sensitity of the biochemcail network as described in the paper'''
    I = np.eye(J.shape[0])
    IJ = np.linalg.inv(I - J)
    for i in range(0, len(IJ)):
        IJ[:, i] = IJ[:, i] / IJ[i, i]
    return(IJ)

def get_data(path, tranpose = False, rem = True):
    '''Read jacobian/model data'''
    data =  pd.read_csv(path, index_col = 0)
    mat = data.values
    if(tranpose):
        mat = mat.T
    if(rem):
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
    
def run_many(path, tranpose = False, rem = True):
    '''Run throw the exisitng models'''
    (data, indexname, df) = get_data(path, tranpose, rem)
    prop = generate_sense_prop(data)
    dist = generate_distance_influence(data)
    neighbors = generate_topology(data)
    return(prop, dist, (neighbors[0], neighbors[1], neighbors[2]))

def compare_mat(mat1, mat2):
    '''Spearman correlation between two list'''
    a = np.linalg.norm(mat1 - mat2)
    b = scipy.stats.spearmanr(np.abs(mat1), np.abs(mat2), axis = None)
    c= scipy.stats.spearmanr((mat1), (mat2), axis = None)

    return((a, b[0], c[0]))


def katz_centrality(mat):
    g = nx.Graph(mat)
    ds = nx.katz_centrality_numpy(g)
    return(ds)
    
    
def example2():
    '''Run a full example'''
    
    '''Get the data and get the sensitvity matricies. '''
    (data, indexname, ind) = get_data(os.path.join("BIOMD0000000313", "jacobian.csv"))
    sensB = senseBio(data)
    (prop, dis, neighbors )= example()
    '''Comptue the correlations'''
    (cor_prob_sd, cor_prop_ud, cor_prob_uu) = (compare_mat(sensB, prop[0] ), compare_mat(sensB, prop[1] + off()),  compare_mat(sensB, prop[2] + off()))
    (cor_dist_ud, cor_dist_uu) = (compare_mat(sensB, dis[0] + off()), compare_mat(sensB, dis[1] + off()))
    (cor_neigh_sd, cor_neigh_ud, cor_neigh_uu) = (compare_mat(sensB, neighbors[0] ), compare_mat(sensB, neighbors[1]) , compare_mat(sensB, neighbors[2] ))    
    
    
    katz1 = katz_centrality(neighbors[1])
    katz1 = np.array(list(katz1.values()))
    katzm1 = (neighbors[1][np.newaxis, :] * katz1 )[0, :, :]
    
    
    katz2 = katz_centrality(neighbors[2])
    katz2 = np.array(list(katz2.values()))
    katzm2 = (neighbors[0][np.newaxis, :] * katz2 )[0, :, :]
    
    (eign0, eign1, eign2) = (np.linalg.eig(neighbors[0])[1], np.linalg.eig(neighbors[1])[1], np.linalg.eig(neighbors[2])[1])
 
    
    '''LAplacian and Eigenvectors'''
    D1 = np.sum(neighbors[2], axis = 1)
        
    D2 = np.sum(neighbors[1], axis = 1)
    L1 = neighbors[2] - D1
    l1 = np.linalg.eig(L1)[1]
    L2 = neighbors[1] - D2
    l2 = np.linalg.eig(L2)[1]
    

    (corrE0, corrE1, corrE2) = (compare_mat(sensB, sens_prop(eign0)), compare_mat(sensB, sens_prop(eign1)), compare_mat(sensB, sens_prop(eign2)))
    (lap1, lap2) = (compare_mat(sensB, sens_prop(katzm1)), compare_mat(sensB, (katzm2)))
   # (lap1, lap2) = (compare_mat(sensB, dist_pop(L1)), compare_mat(sensB, dist_pop(L2)))
    #  D1 = (np.sum(np.abs(), axis = 1))
   
    objects = ('Propogation (d+s)', 'Propogation (d)',
               'Propogation (u)', 'Distance (d)', 'Distance (u)', 
               'Adjacency (d+s)', 'Adjacency (d)', 'Adjacency (u)', 
               'Eigenvectors adj (d+s)', 'Eigvectors adj (d)', 'Eigenvectors adj (u)', \
               'Katz propogation(d)', 'Katz Propogation(u)')
    y_pos = np.arange(len(objects))
    diff_data = [cor_prob_sd, cor_prop_ud, cor_prob_uu, cor_dist_ud, cor_dist_uu, cor_neigh_sd, cor_neigh_ud, cor_neigh_uu, corrE0, corrE1, corrE2, lap1, lap2]
    
    corrs = [x[1] for x in diff_data]
    print(corrs)
    color=iter(cm.rainbow(np.linspace(0,1,len(objects))))
    colors = list(color)
    
    blues = [plt.cm.Blues(300), plt.cm.Blues(200), plt.cm.Blues(100)]
    oranges= [plt.cm.Oranges(200), plt.cm.Oranges(100)]
    greys= [plt.cm.Greys(300), plt.cm.Greys(200), plt.cm.Greys(100)]

    reds= [
            plt.cm.Reds(200), plt.cm.Reds(175), plt.cm.Reds(150),
            plt.cm.Reds(125), plt.cm.Reds(100)
           ]

    colors = blues + oranges + greys + reds
    plt.bar(y_pos, corrs, align='center', alpha=0.5, color = colors)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Correlation with biomodel')
    plt.show()

    corrs2 = [x[2] for x in diff_data]
    plt.bar(y_pos, corrs2, align='center', alpha=0.5, color = colors)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Signed Correlation with biomodel')
    plt.show()
    
    #Calculate the number of equal signed parts. 
    q = lambda inp : np.sum(np.sign(inp) == np.sign(sensB)) / (sensB.shape[0] * sensB.shape[1])
    sign_match = [q(prop[0]), q(prop[1]), q(prop[2]), q(dis[0]), q(dis[1]), 
                  q(neighbors[0]), q(neighbors[1]), q(neighbors[2]), q(eign0), q(eign1), q(eign2), q(l1), q(l2)]
    
    print("SIGN Match")
    print(sign_match)
    plt.bar(y_pos, sign_match, align='center', alpha=0.5, color = colors)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Percent agreement')
    plt.title('Percent agreement in sign')
    plt.show()
    

    #Norm distance
#example2()    
                   
#nx.draw_networkx_nodes(g,pos)
#nx.draw_networkx_edges(g,pos)
#labels = {}
#for i in range(len(ind.columns)):
#    labels[i] = ind.columns[i]
#nx.draw_networkx_labels(g, pos, labels)
#_, v) = np.linalg.eig(neighbors[0])
