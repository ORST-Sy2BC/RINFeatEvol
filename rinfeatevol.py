'''
IMPORTS
'''

import datetime     # for organizing file output
import Bio      # Bio for downloading PDBs, parsing genbank files, printing protein lengths in residues
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio import SeqIO
import pandas as pd
import pypdb
import numpy as np
import os

'''
FUNCTIONS LIST
'''

# Pull coding sequences from a genbank assembly file
def getFeatures(gb_file: str) -> dict:
    '''
    Input: the path to a genbank file for a genome assembly 

    Output: a dictionary containing all CDS feature protein id's as keys and the protein sequence as the value
    '''
    coding_seqs = {}
    for seq_record in SeqIO.parse(gb_file, "genbank"):
        for feature in seq_record.features:
            if feature.type == "CDS":
                coding_seqs[feature.qualifiers['protein_id'][0]] = feature.qualifiers['translation'][0] 
    return coding_seqs

# Query the RCSB PDB (using pypdb) for structures.
def findStrucs(query: str) -> pd.DataFrame:
    '''
    Finds structures matching a RCSB PDB query, and returns the dataframe with their information.
    '''
    search_dict = pypdb.Query(query, query_type="sequence")     # create a dictionary containing search information
    # NOTE: ONLY finds the first 500 for now, to limit download size!
    found = search_dict.search(search_dict)[:500]      # create a list of these PDBs by searching RCSB
    metadata = []           # create a list with the information and the metadata

    for proteins in found:  # for items in # for the items in the list,
        metadata.append(pypdb.describe_pdb(proteins))  # append the dictionary 
        
    return pd.DataFrame(metadata)   # convert, return a Pandas DF

# Sort structures by deposition date.
def sortStrucsByDate(prots: pd.DataFrame) -> pd.DataFrame():
    '''
    Sorts a pandas.DataFrame containing an RCSB PDB query from pypdb by date, and returns the sorted DataFrame.
    '''
    date = [i['deposit_date'] for i in prots.rcsb_accession_info]      # list of deposit dates
    pdbid = [i['id'] for i in prots.entry]
    
    return pd.DataFrame(pdbid,date).sort_index(axis=0)    # return the sorted, reduced dataframe

# Download the set of structures from the query above.
def dlSortedStrucs(prots: pd.DataFrame) -> str:
    '''
    Downloads a set of structures from the above query using the PDB_dl_dir.
    '''
    now = datetime.datetime.now()
    def now_dir_ts():
        '''
        Computes the timestamp for "now", when the query is called
        '''
        now_ts = str(now.year)+"_"+str(now.month)+"_"+str(now.day)+"_"+str(now.hour)+"_"+str(now.minute)+"_"+str(now.second)
        return now_ts

    now = now_dir_ts()      # get the time

    PDB_dl_dir = "ds_"+now  # make the timestamp, save to the class variable

    parser = PDBParser()       # create a parser
    pdbl = PDBList()

    # Download all PDB structures in the previous list if they aren't there
    for pdbid in prots[0]: # index the zeroth col
        pdbl.retrieve_pdb_file(pdb_code=pdbid, file_format='pdb', pdir=PDB_dl_dir)   # Retrieve in PDB format, put in directory 'PDB'

    print('\n#############~DOWNLOAD COMPLETE~#############\n')       # Finished, print "Downloading ... finished!"

    for file in os.scandir(PDB_dl_dir):
        if (file.path.endswith(".ent") and file.is_file()):
            newfn = file.name.replace("pdb","").replace(".ent",".pdb")
            os.rename(file, PDB_dl_dir+"/"+newfn)

    return 

# Walk through the set of structures and obtain a set of lists of all unique proteins in the dataset. Can name them 1, 2, 3 for now.
#    - Hint: using BLAST, unique proteins should be >90% sequence homology with others in the set.
#    - Careful: make sure this algorithm handles possible frameshifts!
#        - start the comparison at the first residues, and as long as thresh residues coincide, the proteins are the same?
def partitionDSbyProtType(path: str) -> list():
    '''
    Takes a downloaded dataset and returns a list of lists, where each inner-list contains protein structures and each outer-list is partitioned by whatever structures are in the files.
    '''

    # for each structure in the file, (BIO pkg)
        # fracture the protein into separate substructures
        # obtain the sequences of each protein as separate items in a holding list
        # if the sequence is < x% homologous to any other existing protein in the outer-list,
            # create a new protein list in the outer-list and add this protein to it
        # else
            # add the current protein to whichever has the greatest sequence homology (if > y%)
    # note: if not working due to frameshifts, try comparing it to the first 20 ++ 10 res / iter 

    return # listname


# Construct the basis matrix for the RIN's
#    - Must be the size of (len(natoms) x len(natoms)) to account for all possible internal interactions!
#    - Number each residue, instead of giving residue names (these might change!)
#    - Be sure to construct the basis sequence carefully using the above functions!!!
#    - Potentially necessary: slice off first ~20 and last ~20 residues to normalize length of the protein before calculating RIN adjacency matrix for some network type. Might make this a parameter of the input function that slices off p % of the front and end of each structure.
#def makeRINcompBasisMat() -> np.array():
    
    
    
#    return 

# Summarize the set of downloaded structures with summary statistics, plot.

# Compute residue interaction network for a protein

# Measure differences between adjacent protein adjacency matrix (i-1) and (i) in the list.
    # Walk through each step, and 

    # Calculate the difference of the adjacency matrix.

# Plot these differences between the residues in the adjacency matrix using cylinders to denote edges.
#   - Draw cylinders (of radius proportional to the absolute value of the difference, scaled appropriately) between nodes with non-zero difference.

# Do the above actions between adjacent proteins, sorted by deposition date.

# Animate this transition.
