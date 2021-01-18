# %%
# imports

import rinfeatevol as rfe

# get the coding sequences from the SARS-CoV-2
id_seq = rfe.getFeatures('sars2.gb')
# %%

# test case if query has zero results
#df = rfe.findStrucs(id_seq['YP_009724393.1'])
#sdf = rfe.sortStrucsByDate(df)

#for key in id_seq:
#    df = rfe.findStrucs(id_seq[key])
#   sdf = rfe.sortStrucsByDate(df)
#    rfe.dlSortedStrucs(sdf) 

#### above functions seem to be working well, however may need more testing. ####
rfe.partitionDSbyProtType('testing/testPDBdirectory', 10.0)
# %%
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
structure = PDBParser().get_structure('6XEY', 'PDB/6XEY.pdb')
ppb=PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())

# %%
from Bio.PDB.StructureAlignment import StructureAlignment

# %%
# extract sequences from the two structures

parser = PDBParser()
structure = parser.get_structure("6XEY", "PDB/6XEY.pdb")[0]['A']

ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())

# %%
# align two structures
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
alignments = aligner.align("ACCGTCAGGAG", "ACGTGACT")
for alignment in sorted(alignments):
    print("Score = %.1f:" % alignment.score)
    print(alignment)
# %%
alignment = pairwise2.align.globalxx("ACCGT", "ACG")
alignment[0][2]

# %%
seq1 = "ACGTGTACTGATGTGTGTGATCACAC"
seq2 = "GTATGATACACACTGATGATGTCATCA"
fracSeqIdentity(seq1,seq2)

# %%
def sameChain(chain1: Bio.PDB.Structure.Structure, chain2: Bio.PDB.Structure.Structure, thresh: float) -> bool:
    '''
    Returns `1` if the two objects passed are the same chain. Returns `0` if they are different chains.
    '''

    def getChainSeq(chain: Bio.PDB.Structure.Structure) -> str:
        '''
        Returns the sequence of a protein chain, in the format of a Bio.PDB.Structure.Structure object.

        ASSUMPTIONS:
        -------------------
        - The chain objects are disjoint (in sequence) from other chain objects.
        - A structure object is not passed (may error-out!)
        '''
        ppb = PPBuilder()
        return str(ppb.build_peptides(chain1)[0].get_sequence())

    from Bio import pairwise2

    def fracSeqIdentity(seq1: str, seq2: str) -> float:
        '''
        Computes the fractional sequence identity between two (non-aligned) sequences. It does this by dividing the score by the greater of the lengths of the two sequences.
        '''
        complength = max(len(seq1), len(seq2))  # compute the max seq length, for normalization
        alignment = pairwise2.align.globalxx(seq1, seq2)   # compute alignment score between the two sequences
        score = alignment[0][2]     # extract the score
        return score/complength     # normalize the alignment score to the larger of the two lengths

    if (fracSeqIdentity(getChainSeq(chain1), getChainSeq(chain2)) > thresh): # if the two chains are above the threshold,
        return True     # return True
    else:               # else,
        return False    # return False

# %%
# extract sequences from a single structure object, then make pairwise comparisons and see which indices are the same on a correlation plot!

items = []

for pp in ppb.build_peptides(structure):
    items.append(str(pp.get_sequence().split("\n")[0]))
items
# %%
# split protein structure into chain structures
import Bio.PDB

def listChains(structure: Bio.PDB.Structure.Structure) -> list:
    '''
    Returns a list of chain objects by splitting the structure object.
    '''
    structures = []
    for chain in structure.get_chains():
        structures.append(chain)
    
    return structures

# %%
parser = PDBParser()
structure = parser.get_structure("6XEY", "PDB/6XEY.pdb")[0]
listChains(structure)

# %% 
# run a pairwise comparison between all chains in the structure

struclist = listChains(structure)

for chain in struclist:
    for chain2 in struclist:
        if (sameChain(chain, chain2, 0.9)):
            print(fracSeqIdentity(getChainSeq(chain), getChainSeq(chain2)))
# %%
seqA = struclist[0]
seqB = struclist[1]

seqA