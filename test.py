# %%
# imports

import rinfeatevol as rfe

# get the coding sequences from the SARS-CoV-2
id_seq = rfe.getFeatures('sars2.gb')
# %%

# test case if query has zero results
df = rfe.findStrucs(id_seq['YP_009724393.1'])
sdf = rfe.sortStrucsByDate(df)

for key in id_seq:
    print(key)
    df = rfe.findStrucs(id_seq[key])
    sdf = rfe.sortStrucsByDate(df)
    print(sdf)
    rfe.dlSortedStrucs(sdf) 


