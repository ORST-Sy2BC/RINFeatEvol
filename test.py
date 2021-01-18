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
# extract sequences from the two structures

parser = PDBParser()
structure = parser.get_structure("6XEY", "PDB/6XEY.pdb")
struclist = listChains(structure)

for chain in struclist:
    for chain2 in struclist:
        if (sameChain(chain, chain2, 0.9)):
            print(rfe.fracSeqIdentity(rfe.getChainSeq(chain), rfe.getChainSeq(chain2)))
