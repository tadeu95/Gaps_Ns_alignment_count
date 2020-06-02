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
        
def cut_alignment(fasta_file,beginning,end):
    """
    Takes in a fasta formatted multiple sequence aligment, a beginning and an end position, then cuts each sequence in the alignment by those positions. Returns a new trimmed fasta file. 
    The first position of the alignment is given by the number 1 in this function.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        record.seq=record.seq[beginning-1:end]
    SeqIO.write(records, "trimmed_alignment.fasta", "fasta") 
   
def first_last_alignment_positions(fasta_file_name,output_file_name):
    '''
    Takes in a multiple sequence alignment in fasta format and the name of the output file and returns a tsv file with three columns: the names of the records in the alignment, the position of the first base pair for that record's sequence and the position of the last base pair for that record's sequence.
    '''
    records = list(SeqIO.parse(fasta_file_name, "fasta"))
    names=[]
    for record in records:
        names.append(record.id)
    sequences=[]
    for record in records:
        sequences.append(record.seq)
    new_list = [str(e) for e in sequences ]
    def find_i(word):
        first = None
        last = None
        for i, letter in enumerate(word):
            if first == None:
                if letter != '-':
                    first = i       
            else:
                if letter != '-':
                        last = i
        return (first, last)
    r = list(map(find_i, new_list))
    first_positions = [i[0]+1 for i in r]
    last_positions = [i[1]+1 for i in r]
    df=pd.DataFrame({'name':names,'first_bp':first_positions,'last_bp':last_positions})
    df.to_csv(output_file_name,sep="\t",index=False)
