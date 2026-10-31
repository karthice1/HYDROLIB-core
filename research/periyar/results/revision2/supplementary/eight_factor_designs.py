import pandas as pd, numpy as np
exec(open('/home/user/ext/ms/an/round4.py').read().split("print('--- subsample")[0].split("for name,r in [('lhs200")[0])
b=pd.read_csv('/home/user/ext/rev2/lhs8b/results.csv').merge(pd.read_csv('/home/user/ext/rev2/lhs8b/samples.csv'),on='run')
both=pd.concat([l8,b],ignore_index=True)
for name,r in [('lhs8 s46',b),('lhs8 pooled 400',both)]:
    a=r.max_flooded_area_km2
    print(name,len(r),'prior',a.min(),a.max(),np.round(np.quantile(a,[.05,.5,.95]),1))
    for t in (0.5,1.0,1.5,2.0):
        print('  tol',t,'both',C(r,acc(r,t)),'HFL',C(r,acc(r,t,radar=False))[0::5])
    print('  radar only',C(r,acc(r,0,hfl=False)))
    k=acc(r,1.0); print('  corr depth-nfp',round(np.corrcoef(r.channel_depth[k],r.n_floodplain[k])[0,1],2),'acc 5-95',np.quantile(a[k],[.05,.95]).round(1))
    print(r.loc[k,['n_channel','q_idamalayar','tide_offset']].agg(['min','max']).round(2).to_string())
