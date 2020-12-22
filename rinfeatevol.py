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

from Bio.SeqUtils import seq1

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
    try:
        search_dict = pypdb.Query(query, query_type="sequence")    # create a dictionary containing search information
        # NOTE: ONLY finds the first 500 for now, to limit download size!
        found = search_dict.search(search_dict)[:500]      # create a list of these PDBs by searching RCSB
        metadata = []           # create a list with the information and the metadata

        for proteins in found:  # for items in # for the items in the list,
            metadata.append(pypdb.describe_pdb(proteins))  # append the dictionary 
        return pd.DataFrame(metadata)   # convert, return a Pandas DF
    except:
        # if no search results are found, return an empty df to be caught in downstream functions.
        print("There were no search results found for the query: " + query)
        return pd.DataFrame()

# Sort structures by deposition date.
def sortStrucsByDate(prots: pd.DataFrame) -> pd.DataFrame:
    '''
    Sorts a pandas.DataFrame containing an RCSB PDB query from pypdb by date, and returns the sorted DataFrame.
    '''

    # if the prots dataframe is empty, no search results were found, simply return the empty df
    if prots.empty:
        return prots
    
    date = [i['deposit_date'] for i in prots.rcsb_accession_info]      # list of deposit dates
    pdbid = [i['id'] for i in prots.entry]
    
    return pd.DataFrame(pdbid,date).sort_index(axis=0)    # return the sorted, reduced dataframe

# Download the set of structures from the query above.
def dlSortedStrucs(prots: pd.DataFrame) -> str:
    '''
    Downloads a set of structures from the above query using the PDB_dl_dir.
    '''
    # check is the prots df is empty, if it is exit the function
    if prots.empty:
        return

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


####################### ABOVE FUNCTIONS IMPLEMENTED! ##################################

# Walk through the set of structures and obtain a set of lists of all unique proteins in the dataset. Can name them 1, 2, 3 for now.
#    - Hint: using BLAST, unique proteins should be >90% sequence homology with others in the set.
#    - Careful: make sure this algorithm handles possible frameshifts!
#        - start the comparison at the first residues, and as long as thresh residues coincide, the proteins are the same?
def partitionDSbyProtType(path: str, α: float) -> list: # α = 0.90
    '''
    Takes a downloaded dataset and returns a list of lists, where each inner-list contains protein structures and each outer-list is partitioned by whatever structures are in the files.

    ARGS:
        path    :   Path to a PDB file with proteins in it.
        α       :   Homology significance level; greater than that value, two proteins are the "same" for structure comparison and will be grouped together; less, they are different/distinct.
    '''
    
    partitioned = []

    def createStruct(filename: str) -> Bio.Entity.Structure.Structure:
        '''
        Creates and returns a structure object from a PDB file, given as the filename/filepath (should work for both!).
        '''
        from Bio.PDB.PDBParser import PDBParser
        parser = PDBParser(PERMISSIVE=1)

        structure_id = filename.replace(old=".pdb", new="")     # replace the file extension with nothing, this is the structure ID!
        return parser.get_structure(structure_id, filename)     # return the structure

    def splitPDBfile(struc: Bio.Entity.Structure.Structure) -> list:
        '''
        Takes a single PDB file and splits the file into a list of structure files based on chain ID. 
        
        TODO: Repeated (identical) structures only exist once.
        '''
        structures = []     # create a list to hold the individual structures
        
        model = struc[0]    # select the first model (generally len(struc) = 1 for most, except for trajectories or NMR models)
        for chain in model: # for each chain in the structure (see the data structure hierarchy of a structure object!)
            structures.append(chain)    # append each chain in the file to the list of structure objects
            # NOTE: we will likely need to refine this. We seek to see multiple levels of interaction, including inter-chain and intra-chain, so doing this in combination with finding inter-domain RINs would be the most comprehensive view!
            # SUGGESTION: Define domains using sequence data, if the above method doesn't work, or if Bio methods won't work.

        return structures

    def strucToSeq(chain: Bio.Entity.Structure.Structure) -> str:
        '''
        Parses a structure object and returns the sequence as a 1-letter AA code.
        '''
        res = list(struc.get_residues())        # residue list from the structure
        seq = ""        # sequence to return later

        for r in res:   # for each residue,
            seq += seq1(r.get_resname())    # append the 3-letter code from each residue name to the sequence string
        return seq

    def computeSeqHomology(chainA: Bio.Entity.Structure.Structure, chainB: Bio.Entity.Structure.Structure) -> float: 
        '''
        Computes the homology between two structures, returns as a float (b/w 0.0 and 1.0)
        '''
        homol = 0.0     # default is zero
        # Existing functions likely work well for this, I just don't know any off the top of my head!

        from Bio import pairwise2

        #pairwise2.format_alignment()

        # compute homology here ...
        
        return homol

    # for each PDB file in the dataset,
    for pdb in path:       # use os module to get filenames???
        # files are already ordered by deposition date, so the list "partitioned" constructed will have a set of lists who all also order the components by deposition date
        
        s = createStruct(pdb)       # create the structure
        chains = splitPDBfile(s)    # split the PDB file into chains
        
        seqs = []
        for chain in chains:
            seqs.append(strucToSeq(chain))    # obtain the sequences of each protein in the structure file

        # compare homology here somewhere?
        # DON'T do pairwise for the whole dataset, that'd be costly. Add one, then compare with the first sequence that was added (assumes the first sequence is representative of the rest of them)




        # if a new protein exists / there isn't one homologous to any of the current chains (perhaps a representative chain, or a random one), (always < α)
            # append a new list to the "partitioned" list to track the evolution of this new structure 
        # else,
            # add the protein to the one with greatest homology

        # note: if not working due to frameshifts, try comparing it to the first 20 ++ 10 res / iter 

    return partitioned


# Construct the basis matrix for the RIN's
#    - Must be the size of (len(natoms) x len(natoms)) to account for all possible internal interactions!
#    - Number each residue, instead of giving residue names (these might change!)
#    - Be sure to construct the basis sequence carefully using the above functions!!!
#    - Potentially necessary: slice off first ~20 and last ~20 residues to normalize length of the protein before calculating RIN adjacency matrix for some network type. Might make this a parameter of the input function that slices off p % of the front and end of each structure.
def makeRINcompBasisMat(seqlist: list()) -> np.array:
    '''
    Takes in a list of presumably identical protein structure types, and constructs a base matrix to enable the sequential comparison of this set of structures. This matrix is of size (natoms x natoms) where the number of atoms is equal to the number in the trimmed, representative length of each protein in the set. Then the sequences are aligned and numbers assigned to the starting sequence number. Thus, a single basis matrix will be populated with RIN information for each identical version in a future function.
    
    Output from this function is constructed by necessarily assuming that the proteins:
    -   Have a similar sequence identity.
    -   Have a similar length.
    -   May have some irregularities especially towards the ends of the sequence. Then the ends may be trimmed a slight amount to enable comparison of the interior residues of the full set of proteins.
    '''
    
    # align the list of sequences
    def alignSeqs(seqlist: list()) -> list():
        # align sequences
        
        return # aligned sequences

    # determine how many residues to trim by creating the start and end position variables as a tuple()
    def detResToTrim(seqlist: list()) -> tuple():
        start = 0
        end = len(seqlist(0))
        
        return (start, end) # define the trimming limits

    # call the above functions to create a basis matrix

    basis_mat = np.array()

    return basis_mat

def constructTrimmedRINmat(minorlist: list()) -> np.array:

    # trim the sequence
    def trimSeq(seq: str, start: int, end: int) -> str:
        trimlist = []

        # trim sequences to be the same length

        return trimlist
    
    # apply trimSeqs() to a single protein structure
    def trimStructure(struc: Bio.PDB.Structure.Structure, start: int, end:int) -> Bio.PDB.Structure.Structure:

        selec = struc
        # from the start to the end positions,
        
        # index the PDB object and sub-select these residues to create a new PDB object

        return selec

    def makeRINmat(struc: Bio.PDB.Structure.Structure, start: int, end: int) -> np.ndarray():
        rinMat = np.ndarray

        # apply above functions, then...
        # apply functions from getcontacts (retrofit, or import from library?)

        return rinMat

    # call above functions
    trimmedRinMat = np.ndarray()

    return trimmedRinMat

def makeRINevolTensor(path: str) -> np.ndarray:
    rinEvolTensor = np.ndarray()

    # for all lists in the "major list",:
        # for all proteins in each "minor list":
            # apply trimStruc sequentially
            # append each new matrix as a numpy ND-array

    return rinEvolTensor

# Summarize the set of downloaded structures with summary statistics, plot.

# Compute residue interaction network for a protein

# Measure differences between adjacent protein adjacency matrix (i-1) and (i) in the list.
    # Walk through each step, and 

    # Calculate the difference of the adjacency matrix.

# Plot these differences between the residues in the adjacency matrix using cylinders to denote edges.
#   - Draw cylinders (of radius proportional to the absolute value of the difference, scaled appropriately) between nodes with non-zero difference.

# Do the above actions between adjacent proteins, sorted by deposition date.

# Animate this transition.
