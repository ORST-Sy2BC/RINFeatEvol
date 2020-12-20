

import pandas as pd
import pypdb

# Query the RCSB PDB (using pypdb) for structures.
def findStrucs(query: string) -> pd.DataFrame:
    '''
    Finds structures matching a RCSB PDB query, and returns the dataframe with their information.
    '''
    prots = pd.DataFrame()  # construct a Pandas DataFrame
    # might use the PD dataframe, ask Jacob to do this
    # include query script from Jacob's implementation mixed w/Antony's
    return prots

# Sort structures by deposition date.
def sortStrucByDate(prots: pd.DataFrame -> pd.DataFrame():
    '''
    Sorts a pandas.DataFrame containing an RCSB PDB query from pypdb by date, and returns the sorted DataFrame.
    '''
    return prots    # return the sorted dataframe

# Download the set of structures from the query above, returns their path in the local directory.

# Walk through the set of structures and obtain a list of all unique proteins. Can name them 1, 2, 3 for now.
#    - Hint: using BLAST, unique proteins should be >90% sequence homology with others in the set.
#    - Potentially necessary: slice off first ~20 and last ~20 residues to normalize length of the protein before calculating RIN adjacency matrix for some network type. Might make this a parameter of the input function that slices off p % of the front and end of each structure.
#    - Careful: make sure this algorithm handles possible frameshifts!

# Construct the basis matrix for the RIN's
#    - Must be the size of (len(natoms) x len(natoms)) to account for all possible internal interactions!
#    - Number each residue, instead of giving residue names (these might change!)
#    - Be sure to construct the basis sequence carefully using the above functions!!!

# Summarize the set of downloaded structures with summary statistics, plot.

# Compute residue interaction network for a protein

# Measure differences between adjacent protein adjacency matrix (i-1) and (i) in the list.
    # Walk through each step, and 

    # Calculate the difference of the adjacency matrix.

# Plot these differences between the residues in the adjacency matrix using cylinders to denote edges.
#   - Draw cylinders (of radius proportional to the absolute value of the difference, scaled appropriately) between nodes with non-zero difference.

# Do the above actions between adjacent proteins, sorted by deposition date.

# Animate this transition.
