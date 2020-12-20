# %%
# imports

import rinfeatevol as rfe

# %%
df = rfe.findStrucs("SARS-CoV-2")
# %%
sdf = rfe.sortStrucsByDate(df)
# %%
rfe.dlSortedStrucs(sdf) # %%

# %%
