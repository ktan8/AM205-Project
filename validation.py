'''Reads suppulemeal data'''

import pandas as pd
import os

path = os.path.join("data", "data_jak2_knockdown.txt")
data = pd.read_csv(path, sep = '\t')
data["RNA1_percent"] = data["sh-JAK2-shRNA1"] / data["shRNA-control"]
data["RNA2_percent"] = data["sh-JAK2-shRNA2"] / data["shRNA-control"]


def label_active(row):
    if(row['RNA1_percent'] >1.2 and row["RNA2_percent"] > 1.2):
        return(1)
    elif(row["RNA1_percent"] < 0.8 and row["RNA2_percent"] < 0.8):
        return(-1)
    return(0)

copy['mat_val'] = copy.apply(label_active, axis = 1)

#Read in the network 
path2 = os.path.join("data", "PDL1_pathway", "PDL1_pathway.matrix.csv")
pathway = pd.read_csv(path2)
pathway = pathway.set_index("Unnamed: 0")



