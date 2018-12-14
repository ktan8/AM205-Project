'''Reads suppulemeal data'''

import pandas as pd
import os
import network_matrices
import numpy as np
import scipy.stats


def get_data():
    '''Read the expression data'''
    path = os.path.join("data", "data_jak2_knockdown.txt")
    data = pd.read_csv(path, sep = '\t')
    data["RNA1_percent"] = data["sh-JAK2-shRNA1"] / data["shRNA-control"]
    data["RNA2_percent"] = data["sh-JAK2-shRNA2"] / data["shRNA-control"]
    
    def label_active(row):
        '''Label whether a given knockout was higher express - huristic is to see 
        if it's above 20% points above/below baseline'''
        if(row['RNA1_percent'] >1.2 and row["RNA2_percent"] > 1.2):
            return(1)
        elif(row["RNA1_percent"] < 0.8 and row["RNA2_percent"] < 0.8):
            return(-1)
        return(0)
    
    data['mat_val'] = data.apply(label_active, axis = 1)
    return(data)
    
def get_pathway():
    '''Read in hand currated network'''
    path2 = os.path.join("data", "PDL1_pathway", "PDL1_pathway.matrix.csv")
    pathway = pd.read_csv(path2)
    pathway = pathway.set_index("Unnamed: 0")
    pathway.index.name = "Gene"
    
    conversion_path = os.path.join("data", "convertModelToExprNames.txt")
    conversions = pd.read_table(conversion_path)
    return(pathway, conversions)

def compute_expressions(conversions, df):
    expression_levels_change = {}
    for elem in conversions: 
        model_name = elem[0]
        targets = elem[1].split(',')
        total_control = 0
        total_expression = 0
        for gene in targets:
            avg_df = df[df["Gene Symbol"] == gene].mean()
            total_control += avg_df['sh-JAK2-shRNA1']
            total_expression += avg_df['shRNA-control']
        change = total_control / total_expression - 1
        if(not np.isnan(change)):
            expression_levels_change[model_name] = change
    return(expression_levels_change)




def get_modes(path):
    (prop, dist, neighbors) = network_matrices.run_many(path2, True)
    (_, gene_names, _) = network_matrices.get_data(path2, True)

# Get the indicies we care about. 
#gene_have = list(copy["Gene Symbol"])
#gene_in_path = list(data.keys())

# Now go and compare the expression with the varios sensititivies in direction. 
def get_val(df, index, gene_names):
    a = (df[df["Gene Symbol"] == gene_names[index]]['mat_val']).sum()
    return(np.sign(a))


def quick_correlation(jak2_index, indicies, mat1, exp):
    print(mat1)
    l1 = np.squeeze(np.asarray(mat1[indicies, jak2_index]))
#    l1 = mat1[jak2_index, indicies].tolist()[0]
    print(l1)
    print(exp)
    return(scipy.stats.spearmanr(np.abs(l1), np.abs(exp)))

def wrapper_corr(jak2_index, indicies, exp):
    return(lambda x : quick_correlation(jak2_index, indicies, x, exp))
    
    
    

data = get_data()
(pathway, conversions) = get_pathway()
elc = compute_expressions(conversions.values, data)

path2 = os.path.join("data", "PDL1_pathway", "PDL1_pathway.matrix.csv")


(prop, dist, neighbors) = network_matrices.run_many(path2, True)
(_, gene_names, _) = network_matrices.get_data(path2, True)

genes_have_all = list(elc.items())
genes_have = [x[0] for x in genes_have_all]
genes_have_val = [x[1] for x in genes_have_all]
indicies = []
jak2_index = -1
for i in range(0, len(gene_names)):
    if(gene_names[i] == "JAK2"):
        jak2_index = i
    if(gene_names[i] in genes_have):
        indicies.append(i)


#corrs = np.zeros(8)
f = wrapper_corr(jak2_index, indicies, genes_have_val)
#

# Need to make sure the diretion is right still. 
corr_data = [f(prop[0].T), f(prop[1].T), f(prop[2].T), f(dist[0].T), f(dist[1].T), 
             f(neighbors[0]), f(neighbors[1]), f(neighbors[2])]
#
#
#f(neighbors[2].T)