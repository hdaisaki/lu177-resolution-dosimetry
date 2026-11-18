"""Regenerate every table of the manuscript and supplement from the results; with --check compare with the published values.
Outputs: results/tables.json
"""
from resdos.paths import ROOT
import os; os.chdir(f'{ROOT}/results')
# Regenerates every table of the short manuscript and supplement from the analysis outputs
# Writes results/tables.json; with --check, compares with results/published/tables_published.json.
import json, os, sys, numpy as np, pandas as pd
J=lambda p: json.load(open(p))
SP=f'{ROOT}/results/spect/'; NE=f'{ROOT}/results/nema/'
order=['Tumour','spleen','Right Kidney','Left Kidney','Liver']
NM={'Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney','Liver':'Liver'}
rec={'Tumour':'f4b','spleen':'f4b','Right Kidney':'f2','Left Kidney':'f2','Liver':'f3'}
sg=lambda x: '0.0' if abs(x)<0.05 else (f'{x:+.1f}').replace('-','−')
fmt=lambda x,d=1: (f'{x:.{d}f}').replace('-','−')
tab={'Tumour':[1916.4519,1809.2385,1864.4781,1713.1579,1396.5609,825.8863],'spleen':[374.4555,428.0064,500.3910,433.7019,314.4332,151.7456],
     'Right Kidney':[494.9367,466.2140,403.4738,326.1331,214.4026,82.6633],'Left Kidney':[489.9914,456.4592,404.1572,325.1385,212.3693,82.3926],
     'Liver':[113.4567,162.4198,137.9349,111.2321,84.2231,46.8034]}
S=pd.read_csv('summary_recommended_models.csv'); R0=pd.read_csv('tac_fit_results_0mm.csv')
X={}
# ---------------- main tables
MASS={'Tumour':'15.6 g','spleen':'122 g','Right Kidney':'115 g','Left Kidney':'137 g','Liver':'1.29 kg'}
GTS={'Tumour':'28.1','spleen':'5.39','Right Kidney':'3.60','Left Kidney':'3.61','Liver':'1.70'}
X['t2']=[[f'{NM[v]} ({MASS[v]})',GTS[v]]+[f'{a:.2f} ± {s:.2f}' for a,s in zip(S[S.VOI==v].sort_values('FWHM_mm').AD_Gy,S[S.VOI==v].sort_values('FWHM_mm').AD_SD_Gy)] for v in order]
sm=J('same_mask_rc.json'); X['t3']=[]
for v in order:
    rc=S[S.VOI==v].sort_values('FWHM_mm').RC_vs_0mm.values
    X['t3'].append([NM[v]]+[f'{x:.3f}' for x in rc[1:]]+[f"{np.max(np.abs(np.array(sm[v]['RC_TIA'])-np.array(sm[v]['RC_AD']))):.3f}"])
fine=J('fine_fwhm.json'); lim=pd.read_csv(f'{ROOT}/results/sphere_model/fwhm_limits.csv'); X['t4']=[]
for v in order:
    th=fine[v]['thr']; X['t4'].append([NM[v]+' (phantom)']+[('> 15' if th[k] is None else fmt(th[k])) for k in ('0.95','0.9','0.8')]+['–'])
f15=lambda x: '> 15' if pd.isna(x) else fmt(x)
for V in [2,4,8,16,32,64,128]:
    c=lim[(lim.voxel_mm>4)&(lim.volume_ml==V)].iloc[0]; f=lim[(lim.voxel_mm<4)&(lim.volume_ml==V)].iloc[0]
    X['t4'].append([f'Sphere {V} mL',f15(c['fwhm_rc_0.95']),f15(c['fwhm_rc_0.90']),f15(c['fwhm_rc_0.80']),f15(f['fwhm_rc_0.90'])])
RF=J(SP+'organ_results_final.json')['organs']; L=J(SP+'local_decomp_final.json'); loo=J(SP+'loo_local.json'); P=J('pvc_summary.json'); X['t5']=[]
for v in order:
    l=L[v]; s=P[v]
    X['t5'].append([NM[v],f"{RF[v]['Q_spect']:.3f}",f"{l['fwhm_local_median']:.1f} ({l['fwhm_range'][0]:.1f}–{l['fwhm_range'][1]:.1f})",f"{l['scale_median']:.2f}",
                    f"{np.mean(np.abs(loo[v]['loo'])):.1f} ({np.mean(np.abs(loo[v]['glob'])):.1f})",
                    f"{sg(s['tiac']['unc'])} / {sg(s['tiac']['gtm'])} / {sg(s['tiac']['rc'])}",
                    f"{sg(s['ad']['unc'])} / {sg(s['ad']['gtm'])} / {sg(s['ad']['rc'])}",
                    f"{sg(s['gtm_ad_rng'][0])} to {sg(s['gtm_ad_rng'][1])}"])
# ---------------- supplementary tables
ho=J('holdout_144.json'); V=J('gt_auc_validation.json'); X['s1']=[]
for v in order:
    r=R0[(R0.VOI==v)&(R0.model==rec[v])].iloc[0]; te=S[S.VOI==v].Teff_72_144_h; e=V[v]['integration']
    oth=[e[m]*100 for m in ['f2','f3','f4b'] if m!=rec[v]]
    X['s1'].append([NM[v],rec[v],f'{te.min():.0f}–{te.max():.0f}' if round(te.min())!=round(te.max()) else f'{te.min():.0f}',f'{S[S.VOI==v].max_rel_residual_pct.max():.1f}',
                    f'{r.AD_CV_pct:.1f}',f'{r.extrapolation_pct:.1f}',sg(ho[v][rec[v]]*100),sg(e[rec[v]]*100),sg(e['model-free']*100),f'{sg(min(oth))} to {sg(max(oth))}',sg(V[v]['dose_conversion']*100)])
X['s2']=[[NM[v]]+[f'{x:.1f}'.replace('-','−') for x in S[S.VOI==v].sort_values('FWHM_mm').RE_vs_GT_pct] for v in order]
rd=J('regression_desc.json'); X['s3']=[]
for v in order:
    q=rd[v]; X['s3'].append([NM[v],f"{q['R2l']:.3f}",f"{q['maxres_l']:.2f}",fmt(q['quad'][0],4),fmt(q['quad'][1],3),fmt(q['quad'][2],1),f"{q['R2q']:.3f}",f"{q['maxres_q']:.2f}"])
ms=J('model_sensitivity.json'); X['s4']=[]
for v in order:
    rcs=np.array([ms[v][m]['RC'] for m in ['f2','f3','f4b']]); res=np.array([ms[v][m]['RE'] for m in ['f2','f3','f4b']])
    X['s4'].append([NM[v]]+[f"{ms[v][m]['RC'][-1]:.3f}" for m in ['f2','f3','f4b']]+[f"{np.max(rcs.max(0)-rcs.min(0)):.3f}"]+[sg(ms[v][m]['RE'][0]) for m in ['f2','f3','f4b']]+[f"{np.max(res.max(0)-res.min(0)):.1f}"])
gc=J('gt_closure.json'); de=J('delineation.json'); ovs=J('osirix_vs_slicer.json')
nom={'Liver':1306.7,'Tumour':15.9,'spleen':124.95,'Right Kidney':111.4,'Left Kidney':133.8}; vol={'Liver':1291.9,'Tumour':15.8,'spleen':121.8,'Right Kidney':110.4,'Left Kidney':133.3}
okey={'Liver':'Liver','Tumour':'Tumour','spleen':'Spleen','Right Kidney':'Right kidney','Left Kidney':'Left kidney'}; X['s5']=[]
for v in order:
    o=ovs[okey[v]]; indep=v in ('Liver','Tumour','spleen'); voi=gc[v]['voi_term']; integ=V[v]['integration'][rec[v]]; conv=V[v]['dose_conversion']
    prod=((1+voi)*(1+integ)*(1+conv)-1)*100
    X['s5'].append([NM[v],f"{vol[v]:.1f} ({sg((vol[v]/nom[v]-1)*100)})",(f"{o['dice']:.3f} / {o['msd']:.2f} / {o['hd95']:.1f}" if indep else 'GT compartment'),
                    sg(voi*100),sg(de[v]['vox_gtvoi']*100),sg(de[v]['delineation']*100),sg(integ*100),sg(conv*100),sg(prod),sg(gc[v]['obs'])])
X['s6']=[[NM[v]]+[sg(x) for x in loo[v]['loo']]+[f"{np.mean(np.abs(loo[v]['loo'])):.1f}",f"{np.mean(np.abs(loo[v]['glob'])):.1f}"] for v in order]
fit=J(SP+'fwhm_fit.json'); U=J(SP+'unified_9.5.json')
X['s7']=[['Total activity in FOV (%)']+[sg((r[2]/r[3]-1)*100) for r in fit['abs_rows']],
         ['Background, measured/model']+[f'{x:.2f}' for x in fit['abs_bg']],
         ['Liver interior, measured/model']+[f'{x:.2f}' for x in fit['abs_lc']]]
for v in order: X['s7'].append([NM[v]+', VOI/true']+[f'{x:.3f}' for x in np.array(U[v]['conc_unc'])/np.array(tab[v])])
X['s8']=[[NM[v],sg(P[v]['gtm_t_rng'][0]),sg(P[v]['tiac']['gtm']),sg(P[v]['gtm_t_rng'][1]),sg(P[v]['gtm1_t']),sg(P[v]['tiac']['rc']),sg(P[v]['rc_ct_t'])] for v in order]
cmp=J(NE+'rc_compare.json'); RN=J(NE+'nema_rc.json'); g=J(NE+'gauss_nema.json')['gauss']
od3d={10:(0.52,0.28),13:(1.21,0.44),17:(2.5,0.57),22:(5.61,0.60),28:(11.57,0.63),37:(26.68,0.69)}
vt={int(r['d']):r for r in cmp['vol_threshold']}; geo={int(r['d']):r for r in RN['results']['CT per sphere | ICF 42.82 cps/MBq']}
peak={37:1.16,28:1.12,22:1.23,17:1.12,13:0.89,10:0.49}; gmap=dict(zip([10,13,17,22,28,37],g)); X['s9']=[]
for d in [10,13,17,22,28,37]:
    X['s9'].append([f'{d}',f'{np.pi/6*d**3/1000:.2f}',f'{od3d[d][0]:.2f} / {od3d[d][1]:.2f}',f"{vt[d]['nvox']} / {vt[d]['rc']:.3f}",f"{geo[d]['nvox']} / {geo[d]['rc_conc']:.3f}",f'{gmap[d]:.3f}',f'{peak[d]:.2f}'])
X['s9'].append(['Fit α (mL) / β','–','2.53 / 0.41','2.46 / 0.40','2.16 / 0.39','–','–'])
d=J('liver_spill_full.json'); ref=d['0.0']['Total']; X['s10']=[]; s0=s15=0
for k,lab in [('Liver','Liver (own activity)'),('Tumour','Tumour'),('Background (incl. lung)','Background (incl. lung)'),('Kidneys','Kidneys'),('Spleen','Spleen')]:
    a=d['0.0'][k]/ref*100; b=d['15.0'][k]/ref*100; s0+=a; s15+=b
    X['s10'].append([lab,f'{a:.2f}',f'{b:.2f}',sg(b-a) if abs(b-a)>=0.05 else '0.00'])
tb=d['15.0']['Total']/ref*100
X['s10']+= [['Unassigned (residual)',f'{100-s0:.2f}',f'{tb-s15:.2f}',sg((tb-s15)-(100-s0))],['Total (GT image)','100.00',f'{tb:.2f}',sg(tb-100)]]
E=J('unified_eval.json')
X['s11']=[[f'{NM[v]}, {lab}']+[sg(x) for x in E[v]['conc_err'][m]] for v in order for m,lab in [('unc','uncorrected'),('gtm','GTM'),('rc','RC')]]
T=J('tia_image_decomp.json')['out']; GT={'Tumour':28.12285,'Spleen':5.382068,'Right kidney':3.598351,'Left kidney':3.599789,'Liver':1.702361}; X['s12']=[]
for v in order:
    o=NM[v]; s=T[o]; bulk=(s['interior (>0.99)']['saam']/GT[o]-1)*100; reg=(s['binary VOI']['saam']/s['interior (>0.99)']['saam']-1)*100
    integ=(s['binary VOI']['fit']/s['binary VOI']['saam']-1)*100; tot=(s['binary VOI']['fit']/GT[o]-1)*100
    gb=(s['binary VOI']['gtm']/s['binary VOI']['fit']-1)*100; gw=(s['occupancy-weighted']['gtm']/s['occupancy-weighted']['fit']-1)*100
    prop=((1+P[v]['tiac']['gtm']/100)/(1+V[v]['integration'][rec[v]])-1)*100
    X['s12'].append([o,sg(bulk),sg(reg),sg(integ),sg(tot),sg(gb),sg(gw),sg(prop),f"{E[v]['self_dose_fraction']:.3f}"])
CS=J('calib_sensitivity.json'); X['s13']=[]
for v in order:
    R=CS['scen'][v]; base=[r for r in R if r['c']==1 and r['s']==1][0]
    rng=lambda k,m: f"{sg(min(r[k][m] for r in R))} to {sg(max(r[k][m] for r in R))}"
    X['s13'].append([NM[v],f"{sg(base['tiac']['gtm'])} ({rng('tiac','gtm')})",f"{sg(base['tiac']['rc'])} ({rng('tiac','rc')})",
                     rng('ad','unc'),f"{sg(base['ad']['gtm'])} ({rng('ad','gtm')})",f"{sg(base['ad']['rc'])} ({rng('ad','rc')})",
                     f"{sum(r['closer_ad']=='GTM' for r in R)} / {len(R)}"])
json.dump(X,open('tables.json','w'),indent=1,ensure_ascii=False)
print('wrote results/tables.json')
if '--check' in sys.argv:
    exp=J(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','results','published','tables_published.json')); bad=0
    for k in exp:
        if exp[k]!=X.get(k):
            bad+=1; print('DIFF',k)
            for a,b in zip(exp[k],X.get(k) or []):
                if a!=b: print('   published',a,'\n   recomputed',b)
    print('all tables identical to the published values' if bad==0 else f'{bad} table(s) differ'); sys.exit(1 if bad else 0)
