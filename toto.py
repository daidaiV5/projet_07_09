import os
import sys
import argparse
import pandas as pd
import math
import numpy as np
from collections import Counter
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

def check_file_path(path):
    """ This function checks that a file exist and return an absolute
        path towards the file if it exists.
        Parameter:
            - path : a string representing a path.
    """
    true_path = os.path.expanduser(path)
    if os.path.isfile(true_path):
        return os.path.abspath(true_path)
    msg = "The path: {} is not a file or does not exist.\n".format(path)
    raise argparse.ArgumentTypeError(msg)


def args_check():
    parser = argparse.ArgumentParser(description= """
     it's a program for analysing of allosteric communication and functional local motions using a structural alphabet
     Install
     ------
     The following software and python packages are required to use this
     	pbxplore 1.4.0
	numpy  1.21.2 
	pandas 1.3.3
	matplotlib  3.4.3  
	seaborn 0.11.2


     Usage
     ------
     for using the toto.py, It will be requires a file fasta who produce by PBxplore.
     
     
     output
     ------
     png:
     
     matrix:matrix de mutual informatique
     """
    ,formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("--fasta", type=check_file_path, help="fasta_produce_by_PBxplore")  
    parser.add_argument("--vitesse_transitions", default=1, required=False, type = int, help="Initial vitesse transitions")
    parameters = parser.parse_args()
    return parameters

def open_fasta(fasta_fichier,vitesse):
    """
    This function for traite the fille fasta, including open the fasta, sorting by order and former the table for the postion de AA and the time of MD simulation
    Parameter:
    fasta_fichier: str,path of fasta
    vitesse: int
    
    output:
    
    """
    dicte={}
    with open(fasta_fichier) as file_one:
        for line in file_one:
            line = line.strip()
            if line.startswith('>'):
                name = line.split(' ', 4 )
                name_model = int(name[3])
                dicte[name_model] = ''
            else:
                dicte[name_model]+=line
    table=pd.DataFrame.from_dict(dicte,orient='index')
    table=table[0].str.split('', expand=True)
    table=table.sort_index()
    table.drop([0,table.shape[1]-1],axis=1,inplace=True)
    liste_index=list(range(1,len(table),vitesse))
    table=table.iloc[liste_index]
    return table

def matrix_transition(transition1,transition2):
    """
    This function create two dictionary of Protein Blocks frequencies at each position(count of each protein/total sequence) and a table of couple of Protein Blocks frequencies at each couple of position(matrix transition for two postion)
    Parameter:
    transition1:position1 from table 
    transition2:postion2 from table
    
    Output:
    test1: dictionnary of frequencies of the transition1
    test2: dictionnary of frequencies of the transition2
    matrix_pourcantage: matrix transition for two postion
    
    """
    sequencelen=len(transition1)
    test1=dict(Counter(transition1))
    test2=dict(Counter(transition2))
    for key,value in test1.items():
        test1[key]=value/sequencelen
    for key,value in test2.items():
        test2[key]=value/sequencelen
    matrix_count=pd.crosstab(pd.Series(transition1,name='position1'),pd.Series(transition2,name='position2'))
    matrix_pourcantage=matrix_count/sequencelen
    return test1,test2,matrix_pourcantage


                                  
def mutual_information(trans,liste1,liste2):
    """
    Compute mutual information between two positions 
    
    Parameters:
    liste1: dictionnary of frequencies of the transition1
    liste2: dictionnary of frequencies of the transition2
    trans: matrix transition for two postion
    
    Output :
    mutual_information: a float value representing mutual information.
    
    """
    MI=0
    for i in list(trans): 
        for j in list(trans.index):
            p_x=liste1[j]
            p_y=liste2[i]
            p_xy=trans.loc[j].at[i]
            if(p_xy != 0):
                MI=MI+p_xy*math.log2(p_xy/(p_x*p_y))
    if(MI>0):
        return MI
    else:
        MI=0
        return MI    

def MI_matrix(table):
    """
    create matrix MI between two position for all postion PB sequence
    
    Parameters:
    table:table of fasta
    
    output:
    matrix MI: pandas.Dataframe containing each mutual information between each postion in a Protein Block sequence. matrix_MI.csv will be stored 
    
    function uesd:
    def matrix_transition
    def mutual_information
    
    """
    MI_dict={}
    for i in range(1,len(table.iloc[0])):
        MI_list=[]
        for j in range(1,len(table.iloc[0])):
            transition1=table[i].tolist()
            transition2=table[j].tolist()
            a,b,c=matrix_transition(transition1,transition2)
            x=mutual_information(c,a,b)
            MI_list.append(x)
        MI_dict[i]=MI_list
    dataframe_MI=pd.DataFrame(MI_dict)
    dataframe_MI.index = list(range(1,len(dataframe_MI)+1))
    for i in list(range(1,len(dataframe_MI))):
        dataframe_MI.iloc[i,i]=0
    dataframe_MI.index.name='position_ver'
    dataframe_MI.columns.name='position_hor'
    dataframe_MI.to_csv("matrix_MI.csv")
    return dataframe_MI    
    
def save_figure(dataframe_MI):
    """ This function plots an heatmap of matrix mutual
        
        Parameters:
            dataframe_MI : pandas.Dataframe containing each mutual information between each postion in a Protein Block sequence.
            
    """
    fig = plt.figure(figsize=(25,25))
    sns.heatmap(dataframe_MI , annot=False, xticklabels=4,yticklabels=4,square = True, linewidths=0.5, cmap="BuGn",)
    plt.tittle=("Matrix Mutal Information")
    plt.savefig("mutual_information_.pdf") 
    
def main():
    parameters = args_check()
    print(parameters) 
    table=open_fasta(parameters.fasta,parameters.vitesse_transitions)
    print('------success for opening fasta----------')
    dataframe_MI=MI_matrix(table)
    save_figure(dataframe_MI)
    print('---------final----------')


if __name__ == "__main__":
    main() 






