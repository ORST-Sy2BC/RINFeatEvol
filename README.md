# RINFeatEvol
Utility to compute the evolution of residue interaction network features over a set of structures along an evolutionary trajectory

# Project title
Evolution of residue interaction network features from single-mutation steps along an evolutionary trajectory

# Algorithm
- Query the RCSB PDB (using pypdb) for SARS-CoV-2 structures.
- Sort structures by deposition date.
- Select all unique proteins in each structure and separate them by type. Collect them as you go.
    - Using BLAST, unique proteins should be >90% sequence homology with others in the set.
- Slice off first ~20 and last ~20 residues to normalize length of the protein before calculating RIN adjacency matrix for some network type.
- Calculate the difference of the adjacency matrix.
- Plot these differences between the residues in the adjacency matrix using cylinders to denote edges.
    - Draw cylinders (of radius proportional to the absolute value of the difference, scaled appropriately) between nodes with non-zero difference.
- Do this between each protein, sorted by deposition date.
- Animate this transition.

# Tool application
- Observe the evolution of a set of COVID proteins, including the spike protein, from structures deposited over the past few months
    - Look for the following things:
        - Proclivity to jump species based on mutation / sequence?
        - Proclivity to become more infectious (mutation D614G?)
