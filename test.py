# %%
# imports

import rinfeatevol as rfe

# get the coding sequences from the SARS-CoV-2
id_seq = rfe.getFeatures('sars2.gb')
# %%

for key in id_seq:
    df = rfe.findStrucs(id_seq[key])
    sdf = rfe.sortStrucsByDate(df)
    rfe.dlSortedStrucs(sdf) 

# %%
