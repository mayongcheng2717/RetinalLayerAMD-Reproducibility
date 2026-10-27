from pathlib import Path
import json, math, shutil, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

BASE = Path('/mnt/agents/output/审稿修订_本地补分析包')
WORK = BASE / 'work_v6'
RR = BASE / 'returned_results'
MR = RR / '09_full54_scdrs_complete/results_for_return/MR'
SC = RR / '09_full54_scdrs_complete/results_for_return/scDRS'
CON = WORK / 'consolidated_results'
FINAL = Path('/mnt/agents/output/Final_Submission_Package')
PNG = FINAL / '03_Figures/PNG'
PDF = FINAL / '03_Figures/PDF'
SVG = FINAL / '03_Figures/SVG'
SOURCE = FINAL / '04_Reproducibility/Source_Data'
for p in [PNG,PDF,SVG,SOURCE]: p.mkdir(parents=True,exist_ok=True)

plt.rcParams.update({
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'path',
    'savefig.facecolor': 'white',
    'savefig.edgecolor': 'white',
})
NAVY='#24557A'; BLUE='#3B7EA1'; TEAL='#2A9D8F'; ORANGE='#E76F51'; GOLD='#E9A23B'
RED='#C94C4C'; GREY='#8B939B'; LIGHT='#E9EEF2'; DARK='#24313B'; GREEN='#4C956C'
LAYER_COLORS={'CSI':'#6D597A','GCC':'#8D6A9F','GCL':'#E76F51','INL':'#2A9D8F','IPL':'#83B5B8','OPL':'#3B7EA1','PS':'#E9A23B','RNFL':'#4C6FFF','RPE':'#C94C4C'}
DISEASE_NAMES={'iamdgc_lateamd':'Late AMD','winkler_earlyamd':'Early AMD','finngen_WET_AMD':'Wet AMD','finngen_DRY_AMD':'Dry AMD','finngen_DM_MACULOPATHY':'Diabetic maculopathy','finngen_H7_RETIVASCOCCLUSION':'Retinal vascular occlusion'}
SHORT_DISEASE={'iamdgc_lateamd':'Late','winkler_earlyamd':'Early','finngen_WET_AMD':'Wet','finngen_DRY_AMD':'Dry','finngen_DM_MACULOPATHY':'DMP','finngen_H7_RETIVASCOCCLUSION':'RVO'}
LAYERS=['RNFL','GCC','GCL','IPL','INL','OPL','PS','RPE','CSI']


def read(p): return pd.read_csv(p)
def export(df,name):
    df.to_csv(SOURCE/f'{name}.csv',index=False)

def axes(fig,rect):
    ax=fig.add_axes(rect)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return ax

def title(ax,text):
    ax.set_title(text,loc='left',fontweight='bold',pad=7,color=DARK)

def label(fig,letter,x,y):
    fig.text(x,y,letter,fontsize=12,fontweight='bold',va='top',ha='left',color='black')

def save(fig,name):
    fig.savefig(PNG/f'{name}.png',dpi=300,facecolor='white')
    fig.savefig(PDF/f'{name}.pdf',dpi=300,facecolor='white')
    fig.savefig(SVG/f'{name}.svg',dpi=300,facecolor='white')
    plt.close(fig)
    print('saved',name)

def clean_axis(ax,xlabel=None,ylabel=None):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.tick_params(length=2.5)

def ptext(p):
    if pd.isna(p): return 'NA'
    if p==0: return '<1e-4'
    if p<0.001: return f'{p:.1e}'
    return f'{p:.3f}'

# Load final result tables.
full54=read(MR/'S3_full54_univariate_complete.csv')
presso=read(MR/'S4_full54_MR_PRESSO_complete.csv')
rev9=read(RR/'03_precision/S4_reverseMR_earlyAMD_precision.csv')
rev_unified=read(CON/'reverseMR_earlyAMD_unified.csv')
mvmr=read(CON/'MVMR_earlyAMD_validated_model.csv')
mvmr_loo=read(CON/'MVMR_earlyAMD_LOO.csv')
mvmr_med=read(CON/'MVMR_earlyAMD_median_sensitivity.csv')
mvmr_infl=read(CON/'MVMR_earlyAMD_influence.csv')
col13=read(RR/'04_smr_coloc/S5a_earlyAMD_H4_gt_0.9_region_audit.csv')
smr42=read(CON/'SMR_legacy42_global_sensitivity.csv')
smrcounts=read(CON/'SMR_global_correction_counts.csv')
med=read(RR/'08_pqtl_grch37_reanalysis/mediation_results.csv').iloc[0]
scmain=read(SC/'main/scdrs_group_summary.csv')
scpar=read(SC/'sensitivity/scdrs_parameter_sensitivity_summary.csv')
sclodo=read(SC/'lodo/scdrs_LODO_summary.csv')
viscorr=read(RR/'05_visium_sensitivity/S6a_Visium_correlation_sensitivity.csv')
vispanel=read(RR/'05_visium_sensitivity/S6d_Visium_16gene_panel_audit.csv')
oldsup=Path('/mnt/agents/output/手稿/补充表_S1-S6.xlsx')
visde=pd.read_excel(oldsup,sheet_name='S6b_Visium_DE')
vissec=pd.read_excel(oldsup,sheet_name='S6c_Visium_section_sensitivity')
overlap=read(CON/'earlyAMD_sample_overlap_bias_scenarios.csv')

# Figure 3: bidirectional and multivariable MR.
fig=plt.figure(figsize=(7.2,8.2),facecolor='white')
ax=axes(fig,[0.25,0.735,0.50,0.185])
keys=[('INL','winkler_earlyamd','IVW'),('OPL','winkler_earlyamd','IVW'),('OPL','winkler_earlyamd','MR-PRESSO corrected'),('GCL','finngen_DRY_AMD','IVW')]
rows=[]
for layer,outcome,method in keys:
    if method=='IVW':
        r=full54[(full54.layer==layer)&(full54.outcome==outcome)].iloc[0]
        rows.append({'layer':layer,'outcome':outcome,'method':method,'beta':r.IVW_beta,'low':r.IVW_ci_low,'high':r.IVW_ci_high,'p':r.IVW_p,'fdr':r.q_full54,'nsnp':r.nsnp})
    else:
        r=presso[(presso.layer==layer)&(presso.outcome==outcome)].iloc[0]
        nout=int(r['n_outliers_bonferroni_0.05'])
        df=int(r.nsnp-nout-1)
        crit=stats.t.ppf(0.975,df)
        rows.append({'layer':layer,'outcome':outcome,'method':method,'beta':r.outlier_corrected_beta,'low':r.outlier_corrected_beta-crit*r.outlier_corrected_se,'high':r.outlier_corrected_beta+crit*r.outlier_corrected_se,'p':r.outlier_corrected_p,'fdr':np.nan,'nsnp':int(r.nsnp-nout)})
f3a=pd.DataFrame(rows);export(f3a,'Figure_3A_forward_MR')
y=np.arange(len(f3a))[::-1]
for i,(yy,r) in enumerate(zip(y,f3a.itertuples())):
    color=LAYER_COLORS[r.layer]
    ax.plot([r.low,r.high],[yy,yy],color=color,lw=1.8,solid_capstyle='round')
    ax.scatter(r.beta,yy,s=48,color=color,edgecolor='white',linewidth=.6,zorder=3,marker='o' if r.method=='IVW' else 's')
    ax.text(1.04,yy,f'{r.beta:.3f} ({r.low:.3f} to {r.high:.3f})',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
ax.set_yticks(y);ax.set_yticklabels([f"{r.layer} to {DISEASE_NAMES[r.outcome]} | {'IVW' if r.method=='IVW' else 'corrected'}" for r in f3a.itertuples()])
ax.axvline(0,color=GREY,lw=.7,ls='--');ax.set_xlim(-.02,.38);clean_axis(ax,'Beta (95% CI)')
title(ax,'Forward MR signals passing full54 FDR')
ax.text(0,-.30,'FDR: INL 0.00133; OPL 0.00941; GCL 0.0374. Corrected OPL estimate excludes rs3138142.',transform=ax.transAxes,fontsize=7,color=GREY)
label(fig,'A',0.055,0.955)

ax=axes(fig,[0.20,0.405,0.55,0.215])
r9=rev9.merge(rev_unified[['disease','layer','fdr_early9_t']],on=['disease','layer'],how='left',validate='one_to_one')
r9['layer']=pd.Categorical(r9.layer,LAYERS,ordered=True);r9=r9.sort_values('layer').reset_index(drop=True)
export(r9,'Figure_3B_reverse_MR')
y=np.arange(len(r9))[::-1]
for yy,r in zip(y,r9.itertuples()):
    color=RED if r.fdr_early9_t<.05 else GREY
    ax.plot([r.IVW_ci95_lower,r.IVW_ci95_upper],[yy,yy],color=color,lw=1.4)
    ax.scatter(r.IVW_beta,yy,s=34,color=color,edgecolor='white',linewidth=.5,zorder=3)
    ax.text(1.03,yy,ptext(r.fdr_early9_t),transform=ax.get_yaxis_transform(),va='center',fontsize=7)
ax.text(1.03,1.04,'FDR',transform=ax.transAxes,fontweight='bold',fontsize=7)
ax.set_yticks(y);ax.set_yticklabels(r9.layer)
ax.axvline(0,color=GREY,lw=.7,ls='--');clean_axis(ax,'Early AMD liability to layer thickness, beta (95% CI)','OCT layer')
title(ax,'Reverse MR across nine layer traits')
label(fig,'B',0.055,0.655)

ax=axes(fig,[0.20,0.075,0.58,0.225])
full=mvmr.copy();full['method']='IVW full model';full=full.rename(columns={'ci95_low':'low','ci95_high':'high'})
med3=mvmr_med[mvmr_med.model=='all_118_instruments'].copy();med3['method']='MVMR-Median';med3=med3.rename(columns={'ci95_low':'low','ci95_high':'high'})
loo3=mvmr_loo[mvmr_loo.omitted_SNP=='rs10737680'].copy();loo3['method']='IVW without rs10737680';loo3=loo3.rename(columns={'ci95_low':'low','ci95_high':'high'})
f3c=pd.concat([full[['layer','method','beta','low','high']],med3[['layer','method','beta','low','high']],loo3[['layer','method','beta','low','high']]],ignore_index=True)
export(f3c,'Figure_3C_MVMR_sensitivity')
methods=['IVW full model','MVMR-Median','IVW without rs10737680'];mcolors=[NAVY,GOLD,ORANGE]
ybase={'INL':2.25,'OPL':1.25,'RPE':.25}
for mi,(method,color) in enumerate(zip(methods,mcolors)):
    sub=f3c[f3c.method==method]
    for r in sub.itertuples():
        yy=ybase[r.layer]+(mi-1)*.18
        ax.plot([r.low,r.high],[yy,yy],color=color,lw=1.4)
        ax.scatter(r.beta,yy,s=35,color=color,edgecolor='white',linewidth=.5,zorder=3)
ax.set_yticks([ybase[x] for x in ['INL','OPL','RPE']]);ax.set_yticklabels(['INL (conditional F=29.5)','OPL (conditional F=39.3)','RPE (conditional F=15.2)'])
ax.axvline(0,color=GREY,lw=.7,ls='--');ax.set_ylim(-.15,2.75);clean_axis(ax,'Multivariable beta (95% CI)')
ax.legend(handles=[Line2D([0],[0],marker='o',color=c,lw=1.4,label=m) for m,c in zip(methods,mcolors)],loc='upper left',frameon=False,ncol=1)
title(ax,'Three-exposure multivariable MR and influence sensitivity')
label(fig,'C',0.055,0.335)
save(fig,'Figure_3')

# Figure 4: colocalization, SMR and TNXB.
fig=plt.figure(figsize=(7.2,10.6),facecolor='white')
ax=axes(fig,[0.19,0.806,0.68,0.145])
cc=col13.copy();cc['row_label']=cc.distinct_coordinate_window_id+' | chr'+cc.chr.astype(str)+' | '+cc.layer
cc['yy']=np.arange(len(cc))[::-1]
export(cc,'Figure_4A_colocalization')
for r in cc.itertuples():
    ax.scatter(r.PP_H4,r.yy,s=46,color=LAYER_COLORS[r.layer],edgecolor='white',linewidth=.6,zorder=3)
ax.set_yticks(cc.yy);ax.set_yticklabels(cc.row_label,fontsize=6.2)
ax.set_xlim(.90,1.005);ax.set_ylim(-.7,len(cc)-.3);clean_axis(ax,'Posterior probability of a shared signal (PP.H4)')
ax.axvline(.9,color=GREY,lw=.7,ls='--')
title(ax,'High-support colocalization: 11 windows, eight intervals')
label(fig,'A',0.045,0.985)

ax=axes(fig,[0.19,0.385,0.60,0.34])
smr42['symbol']=smr42.symbol_legacy42_annotation.fillna(smr42.Gene)
def status(r):
    if bool(r.global_Bonferroni_retained_HEIDI_0_05): return 'Global Bonferroni, HEIDI>0.05'
    if bool(r.global_Bonferroni_retained_HEIDI_0_01): return 'Global Bonferroni, HEIDI>0.01'
    return 'Two-stage catalog only'
smr42['status']=smr42.apply(status,axis=1)
gene_order=smr42.groupby('symbol').size().sort_values(ascending=False).index.tolist()
# Within equal counts use first appearance for a stable biologically interpretable order.
first={g:i for i,g in enumerate(smr42.symbol.drop_duplicates())}
gene_order=sorted(gene_order,key=lambda g:(-sum(smr42.symbol==g),first[g]))
diseases=list(SHORT_DISEASE)
export(smr42,'Figure_4B_SMR_42pairs')
status_color={'Global Bonferroni, HEIDI>0.05':TEAL,'Global Bonferroni, HEIDI>0.01':NAVY,'Two-stage catalog only':'#B9C1C8'}
for r in smr42.itertuples():
    x=diseases.index(r.disease);y=gene_order.index(r.symbol)
    ax.scatter(x,len(gene_order)-1-y,s=52,color=status_color[r.status],edgecolor='white',linewidth=.5)
ax.set_xticks(range(6));ax.set_xticklabels([SHORT_DISEASE[d] for d in diseases])
ax.set_yticks(range(len(gene_order)));ax.set_yticklabels(gene_order[::-1],fontsize=6.4)
ax.set_xlim(-.6,5.6);ax.set_ylim(-.7,len(gene_order)-.3)
ax.set_xlabel('Disease dataset');ax.set_ylabel('Candidate gene')
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
legend_labels=['Bonferroni; HEIDI>0.05','Bonferroni; HEIDI>0.01','Two-stage catalog only']
ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=c,markeredgecolor='white',markersize=6,label=s) for s,c in zip(legend_labels,status_color.values())],loc='upper left',bbox_to_anchor=(1.02,1.0),ncol=1,frameon=False,handletextpad=.3,columnspacing=.8,fontsize=6.2)
title(ax,'SMR catalog and global-correction status')
label(fig,'B',0.045,0.790)

ax=axes(fig,[0.12,0.115,0.36,0.17])
counts=smrcounts.copy()
rule_map={'legacy_within_disease_Bonferroni_HEIDI_0_01':'Two-stage\nHEIDI>0.01','legacy42_with_HEIDI_0_05':'Two-stage\nHEIDI>0.05','global_Bonferroni2354_HEIDI_0_01':'Global Bonf.\nHEIDI>0.01','global_Bonferroni2354_HEIDI_0_05':'Global Bonf.\nHEIDI>0.05','global_BH2354_HEIDI_0_01':'Global BH\nHEIDI>0.01'}
counts=counts[counts.rule.isin(rule_map)].copy();counts['label']=counts.rule.map(rule_map)
counts['order']=counts.rule.map({k:i for i,k in enumerate(rule_map)});counts=counts.sort_values('order',ascending=False)
export(counts,'Figure_4C_SMR_rule_counts')
y=np.arange(len(counts));h=.34
ax.barh(y+h/2,counts.pairs,height=h,color=NAVY,label='Gene-disease pairs')
ax.barh(y-h/2,counts.unique_probes,height=h,color=TEAL,label='Unique genes')
for yy,v in zip(y+h/2,counts.pairs): ax.text(v+1,yy,str(int(v)),va='center',fontsize=6.5)
for yy,v in zip(y-h/2,counts.unique_probes): ax.text(v+1,yy,str(int(v)),va='center',fontsize=6.5)
ax.set_yticks(y);ax.set_yticklabels(counts.label,fontsize=6.5);ax.set_xlim(0,112);clean_axis(ax,'Count')
ax.legend(frameon=False,loc='upper right')
title(ax,'Multiplicity and HEIDI rules')
label(fig,'C',0.045,0.320)

ax=axes(fig,[0.60,0.115,0.30,0.17])
ci_total=stats.t.ppf(.975,31)*med.se_total
f4d=pd.DataFrame([
    ['Step A: GCL to TNXB',med.b_LA,med.b_LA-1.96*med.se_LA,med.b_LA+1.96*med.se_LA,med.p_LA],
    ['Step B: TNXB to dry AMD',med.b_PD,med.b_PD-1.96*med.se_PD,med.b_PD+1.96*med.se_PD,med.p_PD],
    ['Legacy total GCL effect',med.b_total,med.b_total-ci_total,med.b_total+ci_total,med.p_total],
    ['Indirect effect',med.b_indirect,med.b_indirect-1.96*med.se_indirect,med.b_indirect+1.96*med.se_indirect,med.p_indirect],
],columns=['term','beta','low','high','p'])
export(f4d,'Figure_4D_TNXB_mediation')
y=np.arange(4)[::-1]
for yy,r in zip(y,f4d.itertuples()):
    color=ORANGE if r.term=='Indirect effect' else NAVY
    ax.plot([r.low,r.high],[yy,yy],color=color,lw=1.5);ax.scatter(r.beta,yy,s=38,color=color,edgecolor='white',linewidth=.5)
    ax.text(1.03,yy,f'P={r.p:.4f}',transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=6.5,clip_on=False)
d_display_labels=['Step A:\nGCL to TNXB','Step B:\nTNXB to dry AMD','Legacy total\nGCL effect','Indirect effect']
ax.set_yticks(y);ax.set_yticklabels(d_display_labels,fontsize=6.5);ax.axvline(0,color=GREY,lw=.7,ls='--');ax.set_xlim(-.06,.37);ax.set_xticks([0,.1,.2,.3]);clean_axis(ax,'Beta (95% CI)')
title(ax,'Exploratory GCL-TNXB-dry AMD chain')
label(fig,'D',0.52,0.320)
fig.text(.12,.055,'TNXB uses the historical 32-instrument GCL mediation pipeline. The indirect effect is nonsignificant.',fontsize=7,color=GREY)
save(fig,'Figure_4')

# Figure 5: scDRS.
fig=plt.figure(figsize=(7.2,6.8),facecolor='white')
ax=axes(fig,[0.16,0.43,0.68,0.45])
traits=['zekavat_CSI','zekavat_GCC','zekavat_GCL','zekavat_INL','zekavat_IPL','zekavat_OPL','zekavat_PS','zekavat_RNFL','zekavat_RPE','iamdgc_lateamd','winkler_earlyamd']
trait_labels=['CSI','GCC','GCL','INL','IPL','OPL','PS','RNFL','RPE','Late AMD','Early AMD']
cells=['Rod','Cone','BC','AC','HC','NN']
sc=scmain.copy();export(sc,'Figure_5A_scDRS_66')
sc['x']=sc.group.map({c:i for i,c in enumerate(cells)});sc['y']=sc.trait.map({t:i for i,t in enumerate(traits[::-1])})
norm=matplotlib.colors.Normalize(vmin=-5,vmax=5);cmap=plt.get_cmap('RdBu_r')
sizes=20+28*np.clip(-np.log10(sc.assoc_mcp.clip(lower=1e-6)),0,5)
ax.scatter(sc.x,sc.y,s=sizes,c=sc.assoc_mcz,cmap=cmap,norm=norm,edgecolor='white',linewidth=.45)
glob=sc[sc.assoc_fdr_global66<.05]
ax.scatter(glob.x,glob.y,s=sizes[glob.index]+28,facecolor='none',edgecolor='black',linewidth=1.0)
ax.set_xticks(range(6));ax.set_xticklabels(cells)
ax.set_yticks(range(len(traits)));ax.set_yticklabels(trait_labels[::-1])
ax.set_xlim(-.7,5.7);ax.set_ylim(-.8,len(traits)-.2);ax.tick_params(length=0)
for s in ax.spines.values():s.set_visible(False)
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,pad=.03,fraction=.035)
cb.set_label('scDRS mcz');cb.outline.set_visible(False)
ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='white',markeredgecolor='black',markersize=7,label='Black ring: global FDR<0.05')],loc='upper right',bbox_to_anchor=(1.0,1.10),frameon=False)
title(ax,'Primary scDRS associations (66 tests)')
label(fig,'A',0.045,0.955)

ax=axes(fig,[0.30,0.09,0.57,0.25])
g=scmain[scmain.assoc_fdr_global66<.05].copy();g['label']=g.trait.map(dict(zip(traits,trait_labels)))+' | '+g.group
g=g.sort_values('assoc_mcz').reset_index(drop=True);export(g,'Figure_5B_scDRS_global9')
y=np.arange(len(g))
ax.hlines(y,0,g.assoc_mcz,color=TEAL,lw=1.5);ax.scatter(g.assoc_mcz,y,s=42,color=TEAL,edgecolor='white',linewidth=.5)
for yy,r in zip(y,g.itertuples()):ax.text(r.assoc_mcz+.10,yy,f'FDR={r.assoc_fdr_global66:.4f}',va='center',fontsize=6.5)
ax.set_yticks(y);ax.set_yticklabels(g.label,fontsize=6.8);ax.set_xlim(0,5.7);clean_axis(ax,'scDRS mcz')
title(ax,'Nine associations passing global FDR')
label(fig,'B',0.045,0.375)
save(fig,'Figure_5')

# Figure 6: descriptive Visium analyses.
fig=plt.figure(figsize=(7.2,8.2),facecolor='white')
ax=axes(fig,[0.10,0.705,0.36,0.22])
markers=sorted(viscorr.marker_set.unique());palette=plt.get_cmap('tab10').colors
marker_color={m:palette[i%10] for i,m in enumerate(markers)}
export(viscorr,'Figure_6A_Visium_correlations')
for m in markers:
    s=viscorr[viscorr.marker_set==m]
    ax.scatter(s.spearman_r,s.partial_spearman_r,s=28,color=marker_color[m],alpha=.85,edgecolor='white',linewidth=.4,label=m)
ax.plot([-.05,.85],[-.05,.85],color=GREY,lw=.7,ls='--');ax.axhline(0,color=GREY,lw=.6)
ax.set_xlim(-.05,.85);ax.set_ylim(-.07,.85);clean_axis(ax,'Raw Spearman correlation','Technically adjusted partial correlation')
ax.legend(frameon=False,fontsize=5.8,loc='upper left',bbox_to_anchor=(0,-.30),ncol=2,handletextpad=.2,columnspacing=.8)
ax.text(.02,.96,'Empirical FDR<0.05: 0/30 raw; 0/30 partial',transform=ax.transAxes,fontsize=6.5,va='top',color=GREY)
title(ax,'Raw and adjusted spatial correlations')
label(fig,'A',0.035,0.965)

ax=axes(fig,[0.63,0.665,0.31,0.27])
panel=vispanel.copy();panel['symbol']=panel.reported_symbol.fillna(panel.original_panel_symbol)
gene_order=panel.drop_duplicates('candidate_id').symbol.tolist()
sections=panel.drop_duplicates('section').section.tolist()
section_short={s:s.replace('GSM7441725_','').replace('GSM7441726_','').replace('GSM7441727_','').replace('GSM7441728_','').replace('GSM7441729_','') for s in sections}
def statcode(s):
    if s=='evaluable':return 1
    if 'zero variance' in s:return 2
    return 0
panel['code']=panel.evaluation_status.map(statcode);export(panel,'Figure_6B_Visium_panel_audit')
mat=np.zeros((len(gene_order),len(sections)))
for r in panel.itertuples():mat[gene_order.index(r.symbol),sections.index(r.section)]=r.code
from matplotlib.colors import ListedColormap
ax.imshow(mat,aspect='auto',cmap=ListedColormap(['#D9DEE3',TEAL,ORANGE]),vmin=0,vmax=2,interpolation='nearest')
ax.set_xticks(range(5));ax.set_xticklabels([section_short[s] for s in sections],rotation=35,ha='right',fontsize=6)
ax.set_yticks(range(len(gene_order)));ax.set_yticklabels(gene_order,fontsize=5.8)
ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,len(gene_order),1),minor=True)
ax.grid(which='minor',color='white',linewidth=.5);ax.tick_params(which='both',length=0)
for s in ax.spines.values():s.set_visible(False)
ax.legend(handles=[Rectangle((0,0),1,1,color=c,label=l) for c,l in [('#D9DEE3','Not retrieved'),(TEAL,'Evaluable'),(ORANGE,'Zero variance')]],loc='upper center',bbox_to_anchor=(.5,-.25),ncol=3,frameon=False,fontsize=6,handlelength=1,columnspacing=.8)
title(ax,'Sixteen-gene panel audit')
label(fig,'B',0.54,0.965)

ax=axes(fig,[0.16,0.285,0.32,0.28])
vd=visde.sort_values('log2FC_CNV_vs_ctrl').reset_index(drop=True);export(vd,'Figure_6C_Visium_spot_expression')
y=np.arange(len(vd));colors=[ORANGE if p<.05 else GREY for p in vd.fdr]
ax.barh(y,vd.log2FC_CNV_vs_ctrl,color=colors,height=.65)
ax.set_yticks(y);ax.set_yticklabels(vd.gene,fontsize=6);ax.axvline(0,color=DARK,lw=.6);clean_axis(ax,'Spot-level log2 fold difference')
ax.legend(handles=[Rectangle((0,0),1,1,color=ORANGE,label='Parametric FDR<0.05'),Rectangle((0,0),1,1,color=GREY,label='FDR>=0.05')],frameon=False,loc='lower right',fontsize=6)
title(ax,'Descriptive CNV versus control spots')
label(fig,'C',0.035,0.605)

ax=axes(fig,[0.64,0.31,0.29,0.24])
siggenes=vd[vd.fdr<.05].gene.tolist()
vs=vissec[vissec.gene.isin(siggenes)].copy();export(vs,'Figure_6D_Visium_section_ratios')
sec_order=vissec.drop_duplicates('section').section.tolist()
mat=np.full((len(siggenes),len(sec_order)),np.nan)
for r in vs.itertuples():mat[siggenes.index(r.gene),sec_order.index(r.section)]=r.log2_ratio_vs_control_mean
vmax=max(abs(np.nanmin(mat)),abs(np.nanmax(mat)))
im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax,interpolation='nearest')
ax.set_xticks(range(len(sec_order)));ax.set_xticklabels([section_short.get(s,s.split('_')[-1]) for s in sec_order],rotation=35,ha='right',fontsize=6)
ax.set_yticks(range(len(siggenes)));ax.set_yticklabels(siggenes,fontsize=6)
for i in range(len(siggenes)):
    for j in range(len(sec_order)):
        if np.isfinite(mat[i,j]):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=5.4,color='white' if abs(mat[i,j])>.60*vmax else 'black')
ax.tick_params(length=0)
for s in ax.spines.values():s.set_visible(False)
cb=fig.colorbar(im,ax=ax,fraction=.06,pad=.08);cb.set_label('log2 ratio',fontsize=6);cb.outline.set_visible(False);cb.ax.tick_params(labelsize=6)
title(ax,'Section-level descriptive ratios')
label(fig,'D',0.54,0.605)
fig.text(.10,.10,'One CNV donor and one control donor were studied. Disease status and donor identity are completely confounded.',fontsize=7,color=GREY)
save(fig,'Figure_6')

# Supplementary Figure S1: complete forward MR and estimator comparison.
fig=plt.figure(figsize=(7.2,8.4),facecolor='white')
outcomes=['winkler_earlyamd','iamdgc_lateamd','finngen_WET_AMD','finngen_DRY_AMD','finngen_DM_MACULOPATHY','finngen_H7_RETIVASCOCCLUSION']
rects=[[.12,.69,.36,.18],[.58,.69,.36,.18],[.12,.46,.36,.18],[.58,.46,.36,.18],[.12,.23,.36,.18],[.58,.23,.36,.18]]
allrows=[]
for outcome,rect in zip(outcomes,rects):
    ax=axes(fig,rect);s=full54[full54.outcome==outcome].copy();s['layer']=pd.Categorical(s.layer,LAYERS,ordered=True);s=s.sort_values('layer')
    allrows.append(s)
    y=np.arange(len(s))[::-1]
    for yy,r in zip(y,s.itertuples()):
        color=RED if r.q_full54<.05 else NAVY
        ax.plot([r.IVW_ci_low,r.IVW_ci_high],[yy,yy],color=color,lw=.9,alpha=.9);ax.scatter(r.IVW_beta,yy,s=18,color=color,edgecolor='white',linewidth=.3)
    ax.set_yticks(y);ax.set_yticklabels(s.layer,fontsize=5.8);ax.axvline(0,color=GREY,lw=.6,ls='--');ax.tick_params(length=2,labelsize=6)
    ax.set_title(DISEASE_NAMES[outcome],loc='left',fontsize=8,fontweight='bold')
    if outcome in ['finngen_DM_MACULOPATHY','finngen_H7_RETIVASCOCCLUSION']:ax.set_xlabel('IVW beta (95% CI)',fontsize=7)
    ax.set_xlim(-1.1,1.1)
export(pd.concat(allrows),'Supplementary_Figure_S1A_full54')
label(fig,'A',0.035,0.965)
fig.text(.12,.91,'Red points pass full54 FDR. Each panel contains nine layer tests.',fontsize=7,color=GREY)
# Estimator comparison in three small panels.
methods=['IVW','MR-Egger','Weighted median','Simple mode','Weighted mode']
keys=[('INL','winkler_earlyamd'),('OPL','winkler_earlyamd'),('GCL','finngen_DRY_AMD')]
for k,(layer,outcome) in enumerate(keys):
    ax=axes(fig,[.12+k*.30,.045,.25,.105])
    r=full54[(full54.layer==layer)&(full54.outcome==outcome)].iloc[0]
    y=np.arange(5)[::-1]
    for yy,m in zip(y,methods):
        b=r[f'{m}_beta'];lo=r[f'{m}_ci_low'];hi=r[f'{m}_ci_high']
        ax.plot([lo,hi],[yy,yy],color=LAYER_COLORS[layer],lw=1);ax.scatter(b,yy,s=20,color=LAYER_COLORS[layer],edgecolor='white',linewidth=.3)
    ax.set_yticks(y);ax.set_yticklabels(methods,fontsize=5.6);ax.axvline(0,color=GREY,lw=.6,ls='--');ax.set_xlim(-.45,1.02);ax.tick_params(length=2,labelsize=5.8)
    ax.set_title(f'{layer} to {DISEASE_NAMES[outcome]}',fontsize=7,loc='left',fontweight='bold')
    if k==1:ax.set_xlabel('Beta (95% CI)',fontsize=7)
label(fig,'B',0.035,0.180)
save(fig,'Supplementary_Figure_S1')

# Supplementary Figure S2: MVMR diagnostics.
fig=plt.figure(figsize=(7.2,6.3),facecolor='white')
ax=axes(fig,[.12,.64,.28,.25])
f=mvmr[['layer','conditional_F_gencov0']].copy();export(f,'Supplementary_Figure_S2A_conditional_F')
ax.bar(f.layer,f.conditional_F_gencov0,color=[LAYER_COLORS[x] for x in f.layer],width=.6);ax.axhline(10,color=RED,lw=.8,ls='--')
for i,v in enumerate(f.conditional_F_gencov0):ax.text(i,v+1,f'{v:.1f}',ha='center',fontsize=7)
ax.set_ylim(0,46);clean_axis(ax,'Conditional F (gencov=0)')
title(ax,'Instrument strength')
label(fig,'A',.035,.955)
ax=axes(fig,[.57,.47,.37,.42])
inf=mvmr_infl.sort_values('Q_contribution_percent',ascending=False).head(15).sort_values('Q_contribution_percent');export(inf,'Supplementary_Figure_S2B_Q_contributions')
y=np.arange(len(inf));colors=[RED if s=='rs10737680' else NAVY for s in inf.SNP]
ax.barh(y,inf.Q_contribution_percent,color=colors,height=.65);ax.set_yticks(y);ax.set_yticklabels(inf.SNP,fontsize=5.8);clean_axis(ax,'Contribution to joint Q (%)')
title(ax,'Fifteen largest Q contributions')
label(fig,'B',.48,.955)
ax=axes(fig,[.13,.09,.78,.27])
rpe=mvmr_loo[mvmr_loo.layer=='RPE'].sort_values('beta').reset_index(drop=True);export(rpe,'Supplementary_Figure_S2C_RPE_LOO')
x=np.arange(len(rpe))
ax.vlines(x,rpe.ci95_low,rpe.ci95_high,color=LIGHT,lw=.7)
colors=[RED if s=='rs10737680' else NAVY for s in rpe.omitted_SNP]
ax.scatter(x,rpe.beta,s=16,c=colors,edgecolor='white',linewidth=.3,zorder=3)
fullbeta=mvmr[mvmr.layer=='RPE'].beta.iloc[0]
ax.axhline(fullbeta,color=RED,lw=.8,ls='--',label='Full-model estimate');ax.axhline(0,color=GREY,lw=.7)
idx=rpe.index[rpe.omitted_SNP=='rs10737680'][0];val=rpe.loc[idx,'beta']
ax.annotate('Omit rs10737680',xy=(idx,val),xytext=(55,.12),fontsize=7,arrowprops=dict(arrowstyle='->',color=DARK,lw=.7))
ax.set_xlim(-2,len(rpe)+1);ax.set_ylim(-.85,.22);clean_axis(ax,'Omitted instruments ranked by RPE estimate','RPE beta after omission')
ax.legend(frameon=False,loc='lower right')
title(ax,'RPE leave-one-out estimates across 118 instruments')
label(fig,'C',.035,.405)
save(fig,'Supplementary_Figure_S2')

# Supplementary Figure S3: overlap scenarios for INL and OPL.
fig=plt.figure(figsize=(7.2,4.8),facecolor='white')
obsvals=[.1,.2,.5,1.0];layers2=['INL','OPL']
for i,layer in enumerate(layers2):
    for j,obs in enumerate(obsvals):
        ax=axes(fig,[.09+j*.225,.57-i*.42,.17,.25])
        s=overlap[(overlap.layer==layer)&(overlap.absolute_observational_logOR_assumption==obs)]
        for strength,color in [('mean_F',NAVY),('minimum_single_SNP_F',ORANGE)]:
            ss=s[s.strength_parameter==strength].sort_values('relative_overlap_fraction')
            x=ss.relative_overlap_fraction.to_numpy();lo=ss.worst_toward_null_shifted_ci95_low.to_numpy();hi=ss.worst_toward_null_shifted_ci95_high.to_numpy();b=ss.worst_toward_null_shifted_beta.to_numpy()
            ax.fill_between(x,lo,hi,color=color,alpha=.14,linewidth=0);ax.plot(x,b,color=color,lw=1.2,label='Mean F' if strength=='mean_F' else 'Minimum single-SNP F')
            export(ss,f'Supplementary_Figure_S3_{layer}_obs{obs}_{strength}')
        ax.axhline(0,color=GREY,lw=.6,ls='--');ax.set_xlim(0,1);ax.set_ylim(-.02,.27);ax.tick_params(labelsize=6,length=2)
        ax.set_title(f'{layer} | logOR {obs}',fontsize=7,loc='left',fontweight='bold')
        if i==1:ax.set_xlabel('Relative overlap fraction',fontsize=6.5)
        if j==0:ax.set_ylabel('Shifted beta and 95% CI',fontsize=6.5)
label(fig,'A',.025,.965)
fig.legend(handles=[Line2D([0],[0],color=NAVY,lw=1.2,label='Mean F'),Line2D([0],[0],color=ORANGE,lw=1.2,label='Minimum single-SNP F')],loc='upper center',bbox_to_anchor=(.55,.98),ncol=2,frameon=False)
fig.text(.09,.055,'Worst-toward-null scenario bounds. These approximations do not replace non-overlapping replication.',fontsize=7,color=GREY)
save(fig,'Supplementary_Figure_S3')

# Supplementary Figure S4: scDRS parameter and donor sensitivity.
fig=plt.figure(figsize=(7.2,6.7),facecolor='white')
g=scmain[scmain.assoc_fdr_global66<.05].copy();g['label']=g.trait.str.replace('zekavat_','',regex=False).str.replace('iamdgc_lateamd','Late AMD',regex=False).str.replace('winkler_earlyamd','Early AMD',regex=False)+' | '+g.group
labels_order=g.sort_values('assoc_mcz',ascending=False).label.tolist()
designs=['main','sensitivity_w35_10_top500','sensitivity_w35_10_top2000','sensitivity_w10_10_top1000']
rows=[]
for label0 in labels_order:
    base_row=g[g.label==label0].iloc[0]
    rows.append({'label':label0,'design':'main','mcz':base_row.assoc_mcz})
    ss=scpar[(scpar.trait==base_row.trait)&(scpar.group==base_row.group)]
    for r in ss.itertuples():rows.append({'label':label0,'design':r.design,'mcz':r.assoc_mcz})
s4a=pd.DataFrame(rows);export(s4a,'Supplementary_Figure_S4A_scDRS_parameters')
ax=axes(fig,[.25,.56,.59,.32])
for r in s4a.itertuples():
    x=designs.index(r.design);y=labels_order.index(r.label)
    norm4=matplotlib.colors.Normalize(-1,5)
    ax.scatter(x,len(labels_order)-1-y,s=55,c=r.mcz,cmap='RdBu_r',norm=norm4,edgecolor='white',linewidth=.5)
    rgba=plt.get_cmap('RdBu_r')(norm4(r.mcz));luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
    ax.text(x,len(labels_order)-1-y,f'{r.mcz:.2f}',ha='center',va='center',fontsize=5.4,color='white' if luminance<.55 else 'black')
ax.set_xticks(range(4));ax.set_xticklabels(['Main','Top 500','Top 2,000','Window 10/10'])
ax.set_yticks(range(len(labels_order)));ax.set_yticklabels(labels_order[::-1],fontsize=6.5);ax.tick_params(length=0)
for s in ax.spines.values():s.set_visible(False)
title(ax,'Parameter sensitivity for the nine global associations')
label(fig,'A',.035,.955)
ax=axes(fig,[.25,.10,.59,.32])
lrows=[]
for label0 in labels_order:
    base_row=g[g.label==label0].iloc[0]
    ss=sclodo[(sclodo.trait==base_row.trait)&(sclodo.group==base_row.group)].copy();ss['label']=label0;ss['main_mcz']=base_row.assoc_mcz;lrows.append(ss)
s4b=pd.concat(lrows);export(s4b,'Supplementary_Figure_S4B_scDRS_LODO')
for i,label0 in enumerate(labels_order[::-1]):
    ss=s4b[s4b.label==label0];y=i
    vals=ss.assoc_mcz.to_numpy()
    ax.scatter(vals,np.full(len(vals),y),s=22,color=NAVY,alpha=.75,edgecolor='white',linewidth=.3)
    ax.scatter(ss.main_mcz.iloc[0],y,s=42,color=ORANGE,marker='D',edgecolor='white',linewidth=.4,zorder=3)
    n=int((ss.assoc_fdr_descriptive_family<.05).sum())
    ax.text(5.35,y,f'{n}/6 FDR<0.05',va='center',fontsize=6)
ax.set_yticks(range(len(labels_order)));ax.set_yticklabels(labels_order[::-1],fontsize=6.5);ax.set_xlim(1.8,6.3);clean_axis(ax,'scDRS mcz after donor omission')
ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=NAVY,markersize=5,label='Donor omission'),Line2D([0],[0],marker='D',color='none',markerfacecolor=ORANGE,markersize=5,label='Main analysis')],frameon=False,loc='upper left')
title(ax,'Leave-one-donor-out robustness within the same atlas')
label(fig,'B',.035,.465)
save(fig,'Supplementary_Figure_S4')

# Figure 2 is the retained, already verified source graphic. Copy its original PNG and true vector PDF.
shutil.copy2('/mnt/agents/output/图_重排/图2_LDSC遗传结构.png',PNG/'Figure_2.png')
shutil.copy2('/mnt/agents/output/图_重排/图2_LDSC遗传结构.pdf',PDF/'Figure_2.pdf')
subprocess.run(['pdftocairo','-svg',str(PDF/'Figure_2.pdf'),str(SVG/'Figure_2.svg')],check=True)
print('saved Figure_2 from retained vector source')
