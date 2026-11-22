"""Write-up check for E128: the model's per-glomerulus bouton profile (real wiring, APL on) of the lead glomeruli."""
import numpy as np, pandas as pd
from experiments import E128_bouton_groups as E
E119 = E._e119()
U, gl = E.model_boutons('real')
df = pd.DataFrame(U, columns=list(E119.ODORS)); df['g'] = gl
m = df.groupby('g').mean()
tot = m.mean(axis=1)
cvg = m.std(axis=1, ddof=0) / m.mean(axis=1)
nb = df.groupby('g').size()
out = pd.DataFrame({'n_boutons': nb, 'mean': tot.round(2), 'cv_over_odors': cvg.round(2),
                    'top_odor': m.idxmax(axis=1), 'top_share': (m.max(axis=1) / m.clip(lower=0).sum(axis=1)).round(2)})
print(out.sort_values('cv_over_odors', ascending=False).to_string())
for g in ('DM5', 'DM3', 'DM4'):
    print(g, m.loc[g].round(1).sort_values(ascending=False).to_dict())
print('odors', list(E119.ODORS))
