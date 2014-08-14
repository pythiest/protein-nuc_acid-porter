#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pn_parse_lib.py
#  
#  Copyright 2014 <k7n7vi@ivan-UbuVBox>
#  
#  This is a library to basically get a list of NCBI protein GIs and
#  use eUTILS and biopython to automatically retrieve a desired 
#  (upstream) nucleotide sequence from the corresponding CDS in GenBank
#  For now designed to work with "old" RefSeq and GenBank IDs (e.g. YP_)
#  not new WP_ IDs
#  
#  The library contains three main functions:
#  - readGIs reads a list of GIs from file
#  - get_nuc_acc receives a protein GI and retrieves the CDS field
#    position feature (nucleotide seq, plus positions therein)
#    this function is duplicate as SeqIO_get_nuc_acc using the SeqIO
#    functionality to get the same info
#  - parse_nuc_acc gets the CDS field feature and parses it

#imports the Entrez module from the Bio (biopython) library
#Bio is the official name for biopython
#we also set the basic entrez parameters (e.g. email address)
from Bio import Entrez, SeqIO
Entrez.email ="ivan.erill@gmail.com"


#system library for closing the process if on error
import sys

#time library for "sleeping" between NCBI queries
import time


#---------------------------------------------------------------
def read_GIs(file_name = "GI_list.txt"):
    """opens file for reading and gets the GIs into a list, assuming 
    that there is one GI per line.    
    """
    #create file handler for reading file
    try:
        GI_file = open(file_name,"r")
    except (IOError, OSError) as file_open_exception:
        print "The file name provided:", file_name, " does not exist"
        sys.exit()
        
      
    #go through the file, reading GIs into a list
    #we first create an empty GI list
    GIs = []
    
    #we iterate the file as a list of lines
    for line in GI_file:
        #remove \n and empty characters from read line
        #(by default, the iteration reads \n as independent lines)
        #rstrip on a "\n" line will essentially make it an empty string
        line=line.rstrip()
        
        #if line is not empty, then
        if line:
            #if line is numeric
            if (line.isdigit()):
                #append read line to GI list
                GIs.append(line)
                
    #return GIs as list
    return GIs

#---------------------------------------------------------------
def parse_nuc_accession(nuc_acc):
    """parses a nucleotide accession with positions
    basically, two possibilities are contemplated:
    1) forward: ID:pos..pos
    2) reverse: complement(ID:pos..pos) 
    
    returns (using efetch standards):
    - genome accession
    - the strand (1 for forward, 2 for reverse)
    - the sequence start
    - the sequence end
    """
    
    #check for "complement" and assign strand
    #also remove complement from string, and parentheses
    if nuc_acc[0:10]=="complement":
        seq_strand="2"
        nuc_acc=nuc_acc.lstrip("complement(")
        nuc_acc=nuc_acc.rstrip(")")
    else:
        seq_strand="1"
        
    #get genome accession and positions
    gen_acc, gen_pos=nuc_acc.split(":")
    
    #get start and stop positions
    seq_start, seq_stop = gen_pos.split("..")
    
    return gen_acc, seq_strand, seq_start, seq_stop   


#---------------------------------------------------------------
def get_nuc_acc(prot_GI):
    """gets a protein GI and uses eUTILS to navigate NCBI and fetch the
    associated CDS nucleotide accession GI with coordinates
    """
    #get handle to DB object (protein DB)
    prot_handle = Entrez.efetch(db="protein", id=prot_GI, retmode="xml")
    #donwload the handler associated data
    prot_records = Entrez.read(prot_handle)
    prot_handle.close()
    time.sleep(5)  #sleep for 5 seconds
    
    #get the feature list from first element in the records list
    features=prot_records[0]["GBSeq_feature-table"]
    
    #for the features, scan until detecting a GBFeature_key equal to 'CDS'
    #we assume that the CDS feature is unique for a protein record
    for a_feature in features:
        if a_feature["GBFeature_key"]=="CDS":
            feature_qualifiers=a_feature["GBFeature_quals"]
            for f_qual in feature_qualifiers:
                if f_qual["GBQualifier_name"]=="coded_by":
                    genome_loc=f_qual["GBQualifier_value"]
                    return genome_loc
                    
    #return "null" if feature not found
    return None
               

#---------------------------------------------------------------
def SeqIO_get_nuc_acc(prot_GI):
    """gets a protein GI and uses eUTILS to navigate NCBI and fetch the
    associated CDS nucleotide accession GI with coordinates
    Uses SeqIO to parse the protein record and extract the CDS + position
    """
    #get handle to DB object (protein DB) in GB/text format
    prot_handle = Entrez.efetch(db="protein", id=prot_GI, rettype="gb", 
                                retmode="text")
    #donwload the handler associated data with GB parsers
    prot_record = SeqIO.read(prot_handle,"gb")
    prot_handle.close()
    time.sleep(5)  #sleep for 5 seconds

    #the SeqIO parser provides us a neat protein record split into 
    #annotations, features, format, sequence... we use the features list
    #to get to the CDS
    features=prot_record.features
        
    #for the features, scan until detecting a feature type equal to 'CDS'
    #we assume that the CDS feature is unique for a protein record
    for a_feature in features:
        if a_feature.type=="CDS":
            feature_qualifiers=a_feature.qualifiers
            #return with "get" the "coded_by" qualifier 
            #so that it returns none if the qualifier does not exist
            genome_loc=feature_qualifiers.get("coded_by")
            return genome_loc[0]        

    #return "null" if feature not found
    return None
            
#---------------------------------------------------------------
def SeqIO_get_nuc_rec(nuc_acc):
    """gets a list with gen_acc, seq_strand, seq_start and seq_stop
       performs a query to Entrez to retrieve the nucleotide record
       at those positions in FASTA format
       Uses SeqIO read to parse the object
    """
    nuc_handle = Entrez.efetch(db="nuccore", id=nuc_acc[0], 
                                retmode="text", rettype="gb",
                                strand=nuc_acc[1], seq_start=nuc_acc[2],
                                seq_stop=nuc_acc[3])
    nuc_records = SeqIO.read(nuc_handle,"gb")
    nuc_handle.close()
    time.sleep(5)  #sleep for 5 seconds
    #return record in fasta format
    return nuc_records.format("fasta")

#---------------------------------------------------------------
def get_nuc_rec(nuc_acc):
    """gets a list with gen_acc, seq_strand, seq_start and seq_stop
       performs a query to Entrez to retrieve the nucleotide record
       at those positions in FASTA format
    """
    
    #retrieve record with Entrez.read xml
    nuc_handle = Entrez.efetch(db="nuccore", id=nuc_acc[0], 
                                retmode="xml", strand=nuc_acc[1], 
                                seq_start=nuc_acc[2],seq_stop=nuc_acc[3])
    nuc_records = Entrez.read(nuc_handle)
    nuc_handle.close()
    
    time.sleep(5)  #sleep for 5 seconds
    
    #get the GBSeq_definition
    seq_def=nuc_records[0]["GBSeq_definition"]
    #get the sequence
    seq=nuc_records[0]["GBSeq_sequence"]
    #get strand symbol
    strand=""
    if nuc_acc[1]=="2":
        strand="c"
    
    #define fasta header
    fasta_header=">"+nuc_acc[0]+":"+strand+nuc_acc[2]+"-" \
                 +nuc_acc[3]+" "+seq_def
    
    #define fasta field
    fasta_field=fasta_header + "\n" + seq

    #return record in fasta format
    return fasta_field   
    
"""
some explicit code example for efetch/read testing
nuc_handle = Entrez.efetch(db="nuccore", id="NC_010681.1", retmode="text", strand="2", seq_start="4165997",seq_stop="4166596")
nuc_records = Entrez.read(nuc_handle)
    
time.sleep(10)  #sleep for 10 seconds
"""        
    
       

