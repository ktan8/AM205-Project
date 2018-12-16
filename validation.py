'''Reads suppulemeal data'''

import pandas as pd
import os
import network_matrices
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def label_active(row):
    '''Label whether a given knockout was higher express - huristic is to see 
    if it's above 20% points above/below baseline'''
    if(row['RNA1_percent'] >1.2 and row["RNA2_percent"] > 1.2):
        return(1)
    elif(row["RNA1_percent"] < 0.8 and row["RNA2_percent"] < 0.8):
        return(-1)
    return(0)

def get_data_jak2():
    '''Read the expression data'''
    path = os.path.join("data", "data_jak2_knockdown.txt")
    data = pd.read_csv(path, sep = '\t')
    data["RNA1_percent"] = data["sh-JAK2-shRNA1"] / data["shRNA-control"]
    data["RNA2_percent"] = data["sh-JAK2-shRNA2"] / data["shRNA-control"]
    
    data['mat_val'] = data.apply(label_active, axis = 1)
    return(data)
    
def get_data_stat5():
    path3 = os.path.join("data", "data_stat5_knockdown.txt")
    d4 = pd.read_csv(path3, sep = '\t')
    d4['percent_change' ] = d4['data.expr.anno.ctl'] / d4['data.expr.anno.kdSTAT5']
    return(d4)
    
def get_pathway(path):
    '''Read in hand currated network'''
    pathway = pd.read_csv(path)
    pathway = pathway.set_index("Unnamed: 0")
    pathway.index.name = "Gene"
    
    conversion_path = os.path.join("data", "convertModelToExprNames.txt")
    conversions = pd.read_table(conversion_path)
    return(pathway, conversions)

def compute_expressions(conversions, df, gene_col = "Gene Symbol", control_col = "shRNA-control", experiment_col = "sh-JAK2-shRNA1"):
    expression_levels_change = {}
    for elem in conversions: 
        model_name = elem[0]
        targets = elem[1].split(',')
        total_control = 0
        total_expression = 0
        for gene in targets:
            avg_df = df[df[gene_col] == gene].mean()
            total_control += avg_df[experiment_col]
            total_expression += avg_df[control_col]
        change = total_control / total_expression - 1
        if(not np.isnan(change)):
            expression_levels_change[model_name] = change
    return(expression_levels_change)



 
def get_modes(path, tran = True, rem = False,):
    (prop, dist, neighbors) = network_matrices.run_many(path2, tran,rem)
    (_, gene_names, _) = network_matrices.get_data(path2, tran, rem)

# Get the indicies we care about. 
#gene_have = list(copy["Gene Symbol"])
#gene_in_path = list(data.keys())

# Now go and compare the expression with the varios sensititivies in direction. 
def get_val(df, index, gene_names):
    a = (df[df["Gene Symbol"] == gene_names[index]]['mat_val']).sum()
    return(np.sign(a))


def quick_correlation(jak2_index, indicies, mat1, exp):
    print(mat1)
    l1 = np.squeeze(np.asarray(mat1[jak2_index, indicies]))
#    l1 = mat1[jak2_index, indicies].tolist()[0]
    print(l1)
    print(exp)
    print(len(l1))
    print(len(exp))
    return(scipy.stats.spearmanr(np.abs(l1), np.abs(exp)))

def wrapper_corr(jak2_index, indicies, exp):
    return(lambda x : quick_correlation(jak2_index, indicies, x, exp))
    
    
    

def run_example(data, pathway, conversions, elc, knockdown_gene ):
  #  data = get_data()    
    path2 = os.path.join("data", "PDL1_pathway", "PDL1_pathway.matrix.csv")
    
    tran = True
    rem  = True
    
    (prop, dist, neighbors) = network_matrices.run_many(path2, tran, rem)
    (_, gene_names, _) = network_matrices.get_data(path2, tran, rem)
    
    genes_have_all = list(elc.items())
    genes_have = [x[0] for x in genes_have_all]
    genes_have_val = [x[1] for x in genes_have_all]
    indicies = []
    jak2_index = -1
    pdl1_index = pathway.index.get_loc("PD-L1")
    
    for i in range(0, len(gene_names)):
        if(gene_names[i] == knockdown_gene):
            jak2_index = i
        if(gene_names[i] in genes_have):
            indicies.append(i)
    
    
    genes_have_val = np.array(genes_have_val)
    #corrs = np.zeros(8)
    f = wrapper_corr(jak2_index, indicies, genes_have_val)
    #
    
    (eign0, eign1, eign2) = (np.linalg.eig(neighbors[0])[1], np.linalg.eig(neighbors[1])[1], np.linalg.eig(neighbors[2])[1])
    D1 = np.sum(neighbors[2], axis = 1)
            
    D2 = np.sum(neighbors[1], axis = 1)
    L1 = neighbors[2] - D1
    l1 = np.linalg.eig(L1)[1]
    L2 = neighbors[1] - D2
    l2 = np.linalg.eig(L2)[1]
    
    
    # Need to make sure the diretion is right still. 
    corr_data = [f(prop[0]), f(prop[1]), f(prop[2]), f(dist[0]), f(dist[1]), 
                 #f(neighbors[0], f(neighbors[1]), f(neighb)), 
                 ]
    
    blues = [plt.cm.Blues(300), plt.cm.Blues(200), plt.cm.Blues(100)]
    oranges= [plt.cm.Oranges(200), plt.cm.Oranges(100)]
    greys= [plt.cm.Greys(300), plt.cm.Greys(200), plt.cm.Greys(100)]
    
    reds= [
                plt.cm.Reds(200), plt.cm.Reds(175), plt.cm.Reds(150),
                plt.cm.Reds(125), plt.cm.Reds(100)
               ]
    
    colors = blues + oranges + greys + reds
    
    objects = ('Propogation (d+s)', 'Propogation (d)',
                   'Propogation (u)', 'Distance (d)', 
                   'Adjacency (d+s)')
    
    y_pos = np.arange(len(objects))
    
    blues = [plt.cm.Blues(300), plt.cm.Blues(200), plt.cm.Blues(100)]
    oranges= [plt.cm.Oranges(200), plt.cm.Oranges(100)]
    greys= [plt.cm.Greys(300)]
    
    colors = blues + oranges + greys + reds
        
    corr_data = [x[0] for x in corr_data]
    plt.bar(y_pos, corr_data, align='center', alpha=0.5, color = colors)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Correlation')
    plt.title('Strength with experiment')
    plt.show()
    
    def isactive(val):
        if(val > 0.1):
            return(1)
        elif(val < -0.1):
            return(-1)
        return(0)
    
    vfunc = np.vectorize(isactive)
        
    genes_have_val_ind = vfunc(genes_have_val)
    q = lambda inp : np.sum(np.sign(inp[jak2_index, indicies]) == np.sign(genes_have_val_ind)) / (len(indicies))
    sign_match = [q(prop[0]), q(prop[1]), q(prop[2]), q(dist[0]), q(dist[1]), 
                     #q(neighbors[0])
                     ]
        
    
    plt.bar(y_pos, sign_match, align='center', alpha=0.5, color = colors)
    plt.xticks(y_pos, objects, rotation='vertical')
    plt.ylabel('Percent agreement')
    plt.title('Percent agrement in sign with different models')
    plt.show()
    
       
    
    all_data_names = ["Propogation (d+s)", "Progation (d)", "Propogation (u)", "Distance (d)", "Distance (u)", "neighbors (signed)"]
    all_data_models = [prop[0], prop[1], prop[0], dist[0], dist[1], neighbors[0]]
    
    
    
    def write_data(index, all_data, all_data_names, pathway):
        list_arr = []
        for model in all_data:
            list_arr.append((np.squeeze(np.asarray(model[:, index]))))
        df_res = pd.DataFrame(list_arr, index = all_data_names, columns = pathway.columns[1:])
        df_res.index.name = "method"
        return(df_res)
    
    d = write_data(pdl1_index , all_data_models, all_data_names, pathway)
    return(d)
#jak_2 = get_data_jak2()
(pathway, conversions) = get_pathway(os.path.join("data", "PDL1_pathway", "PDL1_pathway.matrix.csv"))
#elc = compute_expressions(conversions.values, jak_2)
#knockdown_gene = "JAK2"
#
#run_example(jak_2, pathway, conversions, elc, knockdown_gene)

stat5 = get_data_stat5()
elc2 = compute_expressions(conversions.values, stat5, 'data.expr.anno.GENE_SYMBOL','data.expr.anno.ctl', 'data.expr.anno.kdSTAT5')
knockdown_gene = "STAT5"
run_example(stat5, pathway, conversions, elc2, knockdown_gene)

#
#
#f(neighbors[2].T)