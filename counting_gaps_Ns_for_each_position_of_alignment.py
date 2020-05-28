# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:47:47 2020

@author: Tadeu
"""
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def Ngap_position_count(alignment_fasta_file):
    """
    Takes in a fasta file containing the result of a multiple sequence alignment, and proceeds to count the number of Ns or gaps for each position of the alignment.
    """
    records = list(SeqIO.parse(alignment_fasta_file, "fasta"))
    list_seqs=[]
    for record in records:
        list_seqs.append(record.seq)
    counted = [v.count('N')+v.count('-')+v.count("n") for v in zip(*list_seqs)]
    return counted
        
def save_as_df(counted,file_name):
    """
    Takes in the list generated with the previous function and converts it to a pandas data frame with two columns: The position of the alignment and the number of gaps or Ns in that position. It saves the data frame as a tsv file.
    """
    df = pd.DataFrame({'counts':counted})
    df['position']=range(1,len(df['counts'])+1)
    df.to_csv(file_name,sep="\t",index=False)
    
def plot_Ngap_counts(counted,image_format):
    """
    Takes in the list generated with the previous function and a file format, then plots the result and saves the figure in a file. The file format is passed as a string and can be png, pdf or svg.
    
    """
    df = pd.DataFrame({'counts':counted})
    df['position']=range(1,len(df['counts'])+1)
    df.plot(kind='scatter',x='position',y='counts',color='blue',s = 0.5)
    plt.axhline(y=0.5,color='red',linestyle="--",linewidth=0.7,alpha=0.6)
    plt.xlabel('Alignment position')
    plt.ylabel('Number of gaps or Ns')
    plt.tight_layout()
    if image_format=="png":
        plt.savefig('counts_plot.png')
    elif image_format=="pdf":
        plt.savefig('counts_plot.pdf')
    elif image_format=="svg":
        plt.savefig('counts_plot.svg')
