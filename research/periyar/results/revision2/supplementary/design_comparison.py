import pandas as pd, numpy as np
R='/home/user/hydrolib-core/research/periyar/results/'
d200=pd.read_csv(R+'v5_lhs200/results.csv').merge(pd.read_csv(R+'v5_lhs200/samples.csv'),on='run')
d400=pd.read_csv(R+'v5_lhs400/results.csv').merge(pd.read_csv(R+'v5_lhs400/samples.csv'),on='run')
m=pd.read_csv(R+'revision2/uncertainty_600_tol0.5/runs.csv')
l8=pd.read_csv(R+'revision2/lhs8_tol1.0/runs.csv')
print(len(d200),len(d400),len(m),len(l8))
def acc(r,t,radar=True,hfl=True):
    a=np.ones(len(r),bool)
    if hfl: a&=(r.hfl_error_m.abs()<=t).values
    if radar: a&=((r.model_flooded_km2>=22)&(r.open_hit_rate>=0.75)).values
    return a
def C(r,mask):
    a=r.max_flooded_area_km2.values; n=mask.sum()
    if n<2: return n,None
    cov=(n-1)/(n+1); lo,hi=np.quantile(a,[(1-cov)/2,(1+cov)/2]); acc_=a[mask]
    return n, round(acc_.min(),1), round(acc_.max(),1), round(lo,1), round(hi,1), round(1-(acc_.max()-acc_.min())/(hi-lo),2)
for name,r in [('lhs200 s42',d200),('lhs400 s43',d400),('pooled600',m),('lhs8 s44',l8)]:
    a=r.max_flooded_area_km2
    print(name,'prior range',a.min(),a.max(),'5-95',np.round(np.quantile(a,[.05,.5,.95]),1),'sd',round(a.std(),1))
    for t in (0.5,1.0,1.5,2.0):
        print('  tol',t,'both',C(r,acc(r,t)),'HFL',C(r,acc(r,t,radar=False))[0::5])
    print('  radar only',C(r,acc(r,0,hfl=False)))
fac=['n_floodplain','channel_depth','canopy_factor','q_rating','q_idukki']
for name,r in [('lhs200 s42',d200),('lhs400 s43',d400)]:
    k=acc(r,1.0); print(name,'accepted param ranges ±1.0:'); print(r.loc[k,fac].agg(['min','max']).round(2).to_string())

print('--- subsample 200 of pooled 600 without replacement, C at ±1.0 m')
rng=np.random.default_rng(1); cs=[]; ns=[]
for i in range(5000):
    s=m.iloc[rng.choice(600,200,replace=False)].reset_index(drop=True); k=acc(s,1.0); r=C(s,k)
    if r[1] is not None: cs.append(r[5]); ns.append(r[0])
cs=np.array(cs); print('C 5-50-95',np.quantile(cs,[.05,.5,.95]).round(2),'P(C>=0.73)',(cs>=0.73).mean().round(3),'n 5-95',np.quantile(ns,[.05,.95]))
# eight-factor accepted range against five-factor like-for-like prior with n=22
a=m.max_flooded_area_km2.values; n=22; cov=(n-1)/(n+1); lo,hi=np.quantile(a,[(1-cov)/2,(1+cov)/2]); print('lhs8 accepted vs 5-factor prior(n=22)',round(1-(210.2-150.0)/(hi-lo),2), lo.round(1),hi.round(1))
print('accepted 5-95 lhs8', np.quantile(l8.max_flooded_area_km2[acc(l8,1.0)],[.05,.95]).round(1), 'pooled', np.quantile(m.max_flooded_area_km2[acc(m,1.0)],[.05,.95]).round(1))
