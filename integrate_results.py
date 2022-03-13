#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:01:50 2022

@author: jason

preceeded by /Users/jason/Dropbox/ELA/prog/AAR_vs_MB_w_map_v2022.py

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import geopandas as gpd
import locale
locale.setlocale(locale.LC_ALL, '')  # Use '' for auto, or force e.g. to 'en_US.UTF-8'
from matplotlib.ticker import FormatStrFormatter

# ----------------------------------- procedure to obtain ensemble
def ensemble(fns,df):
    n_members=18
    ids=np.arange(n_members)+1
    
    dc = {}
    for i in range(len(ids)):
        temp=pd.read_csv(fns[i], delimiter=",")
        dc['df'+str(i)]=temp

    df_concat = pd.concat(dc.values())
    
    temp=dc['df'+str(0)]
    
    by_row_index = df_concat.groupby(df_concat.index)
    
    ensemble_std = by_row_index.std()
    
    ensemble_mean = by_row_index.mean()
    
    ensemble_mean.columns
    
    ensemble_mean['region']=df['region']
    ensemble_mean['name']=df['name']
    ensemble_mean['sector_type']=df['sector_type']
    ensemble_mean['ice_region_type']=df['ice_region_type']
    
    ensemble_std['region']=df['region']
    ensemble_std['name']=df['name']
    ensemble_std['sector_type']=df['sector_type']
    ensemble_std['ice_region_type']=df['ice_region_type']
    
    return ensemble_mean,ensemble_std
# -----------------------------------

# ----------------------------------- procedure to compare AAR and mass balance
# python procedure to compare AAR and mass balance
def AARvsMB(x,y,volume):    
    # slope m and intercept AAR0
    m,AAR0 = np.polyfit(x, y, 1)
    # average AAR
    AAR_mean=np.nanmean(y)
    # imbalance metric alpha
    alpha=AAR_mean/AAR0
    # area to volume scaling exponent gamma
    gamma=1
    # imbalance in Gt
    imbalance=(alpha**gamma-1)*volume
    # sea level equivalent
    SLE=-imbalance/362.

    # correlation R, p-statistic
    R,p=stats.pearsonr(x, y) 
    # x and y to draw regression line
    xx = np.linspace(np.min(x), np.max(x),2)
    yx = m*xx+AAR0
    # confidence statistic
    oneminusp=1-p
    confidence=' = '+str("%.3f"%oneminusp)   
    if round(oneminusp, 3)==1.0:
        confidence='>0.999'
    
    return AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE,xx,yx
# -----------------------------------

ly='x'
wo=0 # 1 if writing out data
batch=1 # set to 1 when generating the ensembple members
prt=0
open_plot=0
plt_map=0
do_effective_AAR0=1 # 1 if adjusting TW AAR0 to a target value for LT sectors
do_plot=0

effective_AAR0_nam=''
if do_effective_AAR0:effective_AAR0_nam='_effective_AAR0'

if_icecaps_name=''
thresh=9

do_SMB=0
plt_name='SMB'
plt_name2='surface mass balance'

if do_SMB==0:
    plt_name='TMB'
    plt_name2='total mass balance'

fn='/Users/jason/Dropbox/1km_grid2/sector_info_v3.csv'
df = pd.read_csv(fn, delimiter=",")
lat = df["lat"]
lon = df["lon"]
area = df["area"]
gl_type = df["gl_type"]
region = df["region"]
name = df["name"]
volume = df["volume"]
name2 = df["name"]

tot_vol=sum(volume)
tot_areal=sum(area)
n_sectors=len(area)-2

valid_flag=[0.]*n_sectors

versionx='20200320'
n_years=20
n_large_enough_sectors=473

SID_index0=1 ; SID_index1=2
ALB_index0=1 ; ALB_index1=2
volume_index0=0 ; volume_index1=1
RCM_index0=1 ; RCM_index1=2

if batch:
    SID_index0=0 ; SID_index1=3
    ALB_index0=0 ; ALB_index1=3
    volume_index0=0 ; volume_index1=1
    RCM_index0=0 ; RCM_index1=2

fns=[]
AAR1fns=[]
AAR0fns=[]
MBfns=[]

# ------------------------------------------------ loop over SID uncertainty
devname=["lo","mid","hi","const"]

SID_devname=["0.9","1.0","1.1","const"]

for SID_dev_index in range(SID_index0,SID_index1):    
    # ------------------------------------------------ loop over volume treatments
    volume_name=['unscaled','scaled']
    for volume_index in range(volume_index0,volume_index1):
        volume = df["volume"]
        if volume_index:
            v=np.where(area<89000)
            volume[v[0]]=area[v[0]]*1.6233623010
        namex = df["name"]
        valid_flag=[0.]*n_sectors
        # ------------------------------------------------ loop over RCMs
        RCM=['MAR','RACMO']
        for mm in range(RCM_index0,RCM_index1):
        # ------------------------------------------------ loop over albedo uncertainty
            for ALB_dev_index in range(ALB_index0,ALB_index1):
                fn='/Users/jason/Dropbox/ELA/stats/TSL/'+\
                            namex[0]+'_ELA_v'+versionx+'_'+devname[ALB_dev_index]+'.csv'
                TMB_annual_tot = pd.read_csv(fn, delimiter=",")
                # https://www.dropbox.com/sh/ckn2cnf7b67if61/AABngYeaICreT8j6RhFOqHwxa?dl=0
                TMB_annual_tot=TMB_annual_tot.drop(columns=TMB_annual_tot.columns[1:])
                TMB_annual_tot["TMB"]=0.
                
                lost_count=0
                valid_vol_sum=0.
                valid_area_sum=0.

                # if wo:
                out_fn='/Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/'+\
                        'ALB'+devname[ALB_dev_index]+'_'+RCM[mm]+'_'+\
                            volume_name[volume_index]+'_SID'+SID_devname[SID_dev_index]+\
                                if_icecaps_name
                
                fns.append(out_fn+'.csv')
                AAR1fns.append(out_fn+'.AARs1.csv')
                AAR0fns.append(out_fn+'.AARs0.csv')
                MBfns.append(out_fn+'.MBs.csv')

                # AARs = np.fromfile(out_fn+'2x473x20.AARs.npy', dtype='float32')
                # AARs=AARs.reshape(2,n_large_enough_sectors,n_years)
                # MBs = np.fromfile(out_fn+'1x473x20.MBs.npy', dtype='float32')
                # MBs=MBs.reshape(n_large_enough_sectors,n_years)

                df=pd.read_csv(out_fn+'.csv')
                # print(out_fn)
                # k=2
                # plt.scatter(MBs[k,:],AARs[0,k,:])
                # plt.scatter(MBs[k,:],AARs[1,k,:])
                # plt.title(df.name[k])

##%% concat ensemble members, sector scale AAR and MB

AAR,AARstd=ensemble(AAR1fns,df)
MB,MBstd=ensemble(MBfns,df)
# AAR,AARstd=ensemble(AAR0fns,df)

# #%% diagnostic
# k=2
# plt.title(df.name[k])
# plt.scatter(MB.iloc[k,v],AAR.iloc[k,v])

#%% regional results

from matplotlib.ticker import StrMethodFormatter
formatx='{x:,.2f}'
fn='/Users/jason/Dropbox/1km_grid2/sector_info_v3.csv'
df = pd.read_csv(fn, delimiter=",")

# plt.rcParams['font.sans-serif'] = ['Georgia']
th=2 # line thickness
font_size=28
plt.rcParams["font.size"] = font_size
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = False
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'

regions=['SW','CW','NW','NO','NE','CE','SE','all'] ; n_regions=len(regions)
regions2=['Southwest (SW)','Central West (CW)','Northwest (NW)','North (NO)','Northeast (NE)','Central East (CE)','Southeast (SE)','All']

regions=['SW'] ; n_regions=len(regions)
regions2=['Southwest (SW)']

# regions=['CW'] ; n_regions=len(regions)
# regions2=['Central West (CW)']

regions=['NW'] ; n_regions=len(regions)
regions2=['Northwest (NW)']

regions=['NO'] ; n_regions=len(regions)
regions2=['North (NO']

regions=['all'] ; n_regions=len(regions)
regions2=['all ice']

# regions=['NO'] ; n_regions=len(regions)
# regions=['SW','CW','NW','NO','NE','CE','SE'] ; n_regions=len(regions)
# regions=['SW','CW'] ; n_regions=len(regions)

sector_scenario=0
n_years=20
v=np.arange(n_years)+1 # valid cases

# colors=['r','m','k','b','darkorange','c','grey']

tot_SLR=0
cc=0
do_plot=1

SLR_cum=[]
vol_cum=[]
volsum2=0.
SLRcum2=0.
AAR0_collection=[]
area_collection=[]
name_collection=[]
SLE_collection=[]

for region_index,region in enumerate(regions):
    integrator_sum=0.
    if do_plot:
        plt.close()
        fig, ax = plt.subplots(figsize=(10,10))
    regionx=df.region.values
    print(region)
    AAR_integrated=np.zeros(n_years)
    AAR_err_integrated=np.zeros(n_years)
    MB_integrated=np.zeros(n_years)
    MB_err_integrated=np.zeros(n_years)
    vol_sum=0.
    if region!='all':
        N=1
        tot_regional_area=np.sum(df.area[df.region==region])#[0:N])
        tot_regional_vol=np.sum(df.volume[df.region==region])#)[0:N])
        regionx=df.region.values
    else:
        tot_regional_area=np.sum(df.area)
        tot_regional_vol=np.sum(df.volume)        
        regionx[:]='all'
    for k in range(n_large_enough_sectors):
        if regionx[k]==region:
        # if k>470:
        # if k==4
        # if name[k]=='OSTENFELD_GLETSCHER':
        # # if region[k]=='NO':
        # if df.name[k]=='STORSTROMMEN': 
        # if df.name[k]=='AB_DRACHMANN_GLETSCHER_L_BISTRUP_BRAE': 
        # if df.name[k]=='JAKOBSHAVN_ISBRAE':
        # if df.name[k]=='SAQQAP-MAJORQAQ-SOUTHTERRUSSEL_SOUTHQUARUSSEL':
        # if df.region[k]=='NW':
        # if ((df.region[k]=='SW')&(df.sector_type[k]=='LT')):
            if sector_scenario==0:
                # if regionx[k]==region:
                if regionx[k]!='null':
                    scenario_name='all_TMB'
                    integrator=df.area[k]/tot_regional_area
                    integrator_sum+=integrator
                    temp=AAR.iloc[k,v].values.astype(float)
                    count_of_invalid_AARs=sum(np.isnan(temp))
                    # print(k,count_of_invalid_AARs,temp)
                    if count_of_invalid_AARs>0:
                        # print(df.name[k])
                        # print(k,'before')
                        # print(temp)
                        temp[np.isnan(temp)]=np.nanmean(temp)
                        # print('after')
                        # print(temp)               
                    if count_of_invalid_AARs==n_years:
                        temp[:]=0.5
                        # print(temp)                        

                    if count_of_invalid_AARs<20:
                        x=MB.iloc[k,v].values.astype(float)
                        y=np.array(temp)
                        # print(k,x,y)
                        
                        # AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE,xx,yx=AARvsMB(x,y,df.volume[k])
                        
                        # AAR0_collection.append(AAR0)
                        # area_collection.append(df.area[k])
                        # name_collection.append(df.name[k])
                        # SLE_collection.append(SLE)
                        # volsum2+=df.volume[k]
                        # SLRcum2+=SLE

                        # SLR_cum.append(SLRcum2)
                        # vol_cum.append(volsum2)
                        # print(df.name[k],SLE)
                        # print(temp)

                        cc+=1

                    temp_err=AARstd.iloc[k,v].values.astype(float)
                    count_of_invalid_AARs=sum(np.isnan(temp_err))
                    # print(k,count_of_invalid_AARs,temp)
                    if count_of_invalid_AARs>0:
                        # print()
                        # print(k,'before')
                        # print(temp)
                        temp_err[np.isnan(temp_err)]=np.nanmean(temp_err)
                        # print('after')
                        # print(temp)            

                    if count_of_invalid_AARs==n_years:
                        temp_err[:]=0.025        

                    # integrator=1.
                    AAR_integrated+=np.array(integrator*temp)
                    AAR_err_integrated+=np.array(integrator*temp_err)
                    MB_integrated+=MB.iloc[k,v].values.astype(float)
                    MB_err_integrated+=MBstd.iloc[k,v].values.astype(float)
                    vol_sum+=df.volume[k]

                    # if cc>N-1:
                    #     break                    

                    # y=temp ; AAR_integrated=temp

                    # x=MB_integrated ; y=AAR_integrated
                    # m,AAR0 = np.polyfit(x, y, 1)
                    # R,p=stats.pearsonr(x, y) 

                    # alpha=np.nanmean(y)/AAR0
                    # gamma=1
                    # imbalance_mean=(alpha**gamma-1)*df.volume[k]
                    # SLR=-imbalance_mean/362.
                    # tot_SLR+=SLR
                    # AAR_mean=np.nanmean(y)

                    print("area %.0f" %df.area[k],'count',cc,"AAR0 %.4f" %AAR0,"Area fraciton %.4f" %integrator,"%.4f" %integrator_sum,df.name[k])#,"%.3f" %SLE,"%.3f" %tot_SLE)
                    # print(vol_sum)

    plt.plot(SLR_cum,'-o')
    # print(np.sum(SLR_cum),cc)
    # plt.plot(vol_cum,'-o')
    # print(np.sum(vol_cum))
    do_plot=1
    
    if do_plot:
    
        AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE1,xx,yx=AARvsMB(MB_integrated+MB_err_integrated,AAR_integrated+AAR_err_integrated,vol_sum)
        print('SLE +',SLE1)
        AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE0,xx,yx=AARvsMB(MB_integrated,AAR_integrated-AAR_err_integrated,vol_sum)
        print('SLE -',SLE0)
        mm_err=abs(SLE1-SLE0)
        # AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE,xx,yx=AARvsMB(MB_integrated+MB_err_integrated,AAR_err_integrated-AAR_err_integrated)
        # print('SLE +-',SLE)
        # AAR0,alpha,R,confidence,imbalance,SLE,xx,yx=AARvsMB(MB_integrated-MB_err_integrated,AAR_err_integrated+AAR_err_integrated)
        # print('SLE -+',SLE)
    
        x=MB_integrated ; y=AAR_integrated
        AAR0,AAR_mean,alpha,R,confidence,imbalance,SLE,xx,yx=AARvsMB(x,y,vol_sum)
    
        # txt=" AAR0: %.2f" %AAR0+"\nR: %.2f" %R+"\nSLR: %.1f" %tot_SLR+" mm"
        # txt=" AAR0: %.2f" %AAR0+"\nSLR: %.0f" %tot_SLR+" mm"
    
        ax.scatter(x, y, s=1200, facecolors='k', edgecolors='k',linewidth=th,zorder=9)
    
        plt.errorbar(x, y, yerr=[AAR_err_integrated, AAR_err_integrated],
                     xerr=MB_err_integrated,fmt='.',c='w', ecolor='purple', capthick=30,capsize=10, 
                    elinewidth=th,
                    markeredgewidth=th)
    
        xos=np.std(x)/2
        xos_text=xos/2
        yos_text=0.004
        yos=np.std(y)
        year=np.arange(2000,2021).astype(str)
        for i in range(0,20):
            plt.text(x[i], y[i]-yos_text/10, str(year[i])[2:4], \
            horizontalalignment='center',verticalalignment='center', \
            color='w',zorder=10,fontsize=28)
    
        ax.plot(xx, yx, '--',color='grey',linewidth=th*2)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.set_ylabel('Accumulation Area Ratio (AAR)')
        ax.set_xlabel('mass balance, Gt y$^{-1}$')
        ax.set_ylim(np.min(y)-np.max(AAR_err_integrated),np.max(y)+yos/2)
    
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter(formatx))
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
    
        yx = [AAR0,AAR0]
        xx=np.linspace(0, 0+xos/2,2)
    
        plt.plot(xx, yx, color='b',linewidth=th*2,zorder=20)
    
        if region=='CE':yos_text*=4
        if region=='NE':yos_text*=2
        if region=='NW':yos_text*=2
        if region=='NO':yos_text*=2
    
        plt.text(xx[1]+xos_text,AAR0-yos_text,'AAR$_0$\n'+str("%8.3f"%AAR0).lstrip(),color='b',zorder=12)
        
        yx = [AAR_mean,AAR_mean]
        plt.plot(xx, yx, color='r',linewidth=th*2,zorder=20)
        #                    plt.text(np.min(x)+xos,yx[0]+yx[0]*0.005,'AAR$_0$ = '+str("%8.3f"%AAR0).lstrip(),color='b')
    
        # nam='AAR$_{2000-2019}$'
        txt='$\overline{AAR}$'
        plt.text(xx[1]+xos_text,AAR_mean-yos_text,txt+'\n'+str("%8.3f"%AAR_mean).lstrip(),color='r')
        
        props = dict(boxstyle='round', facecolor='w',edgecolor='w',linewidth=th,alpha=0.7)
    
        # fig = plt.gcf()
        # ax = plt.add_subplot(111)
    
    
        cc=1    ;   xx0=0.03    ;   yy0=0.95    ;   dy=0.06
        ax.text(xx0, yy0+cc*dy,regions2[region_index],color='k',
            transform=ax.transAxes,#, fontsize=9
            verticalalignment='top', bbox=props,fontsize=font_size*1.4)    ;   cc+=1
    
        # txt='AAR$_{2000-2019}$/AAR$_{0}$ = '+str("%8.3f"%alpha_mean).lstrip()+'±'+str("%8.3f"%mean_y_err).lstrip()+\
        #     '\nSLR = '+str("%8.0f"%mm_SLE_mean).lstrip()+'±'+str("%8.0f"%mm_err).lstrip()+' mm'+\
        #     '\nCorrelation = '+str("%8.3f"%R).lstrip()+', 1-p '+confidence
    
        mult=1
        cc=1    ;   xx0=0.45   ;   yy0=0.16    ;   dy=0.06
    
        if region=='NO':xx0=0.6
        if region=='NE':xx0=0.65
        txt='SLR = '+str("%8.0f"%SLE).lstrip()+'±'+str("%8.0f"%mm_err).lstrip()+' mm'+\
            '\ncorrelation = '+str("%8.3f"%R).lstrip()+'\n1-p '+confidence
        props = dict(boxstyle='round', facecolor='w',edgecolor='grey',linewidth=th,alpha=1)
    
        ax.text(xx0, yy0+cc*dy,txt,color='k',
            transform=ax.transAxes,#, fontsize=9
            verticalalignment='top', bbox=props,fontsize=font_size*mult)    ;   cc+=1
    
        ymin=0.5
        yos=0.002
        if np.min(y)-yos < ymin:ymin=0.
        
        xx = [0,0]
        #                yx = [np.min(y),AAR0]
        yx = [ymin,AAR0]
        plt.plot(xx, yx, '--', color='g',linewidth=th*2)
        
        propsk = dict(boxstyle='round', fc='k',facecolor='k',edgecolor='grey',alpha=1)
    
        du=0
        if du:
            cc=0    ;   xx0=0.72    ;   yy0=0.3    ;   dy=0.06
            plt.text(xx0, yy0+cc*dy,'2 digit year',color='w',
                transform=ax.transAxes,#, fontsize=9
                verticalalignment='top', bbox=propsk)       
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        fig_path='/Users/jason/Dropbox/ELA/Figs/AAR_vs_MB_regional_integration/'
        by_area_name='_by_area' ; scenario_name='all_TMB'
        figname=fig_path+region+'_'+scenario_name+by_area_name+if_icecaps_name
        ly='x'
        if ly=='p':
            plt.savefig(figname+'.png', bbox_inches='tight', dpi=300)
            # plt.savefig(figname+'.eps', bbox_inches='tight')
print('tot_SLR',tot_SLR)
#%% concat ensemble members, regional 
n_members=18
ids=np.arange(n_members)+1

dc = {}
for i in range(len(ids)):
    temp=pd.read_csv(fns[i], delimiter=",")
    varx=['AAR_mean_truncated','AAR0_truncated', 'alpha_mean_truncated', 'mmSLE_mean_truncated','AAR_mean', 'AAR0', 'alpha_mean', 'mmSLE_mean', 'AAR_2012','alpha_2012', 'mmSLE_2012', 'AAR_2018', 'alpha_2018', 'mmSLE_2018']
    for var in varx:
        temp[var][temp['valid_flag'] == 0]=np.nan
    dc['df'+str(i)]=temp
    # print(dc['df'+str(i)])
df_concat = pd.concat(dc.values())

temp=dc['df'+str(0)]

#%%
plt.plot(temp.mmSLE_2012)

#%% ontain ensemble
df=temp
by_row_index = df_concat.groupby(df_concat.index)

df_ensemble_std = by_row_index.std()

df_ensemble_mean = by_row_index.mean()

df_ensemble_mean.columns
df.columns

df_ensemble_mean['region']=df['region']
df_ensemble_mean['name']=df['name']
df_ensemble_mean['sector_type']=df['sector_type']
df_ensemble_mean['ice_region_type']=df['ice_region_type']

df_ensemble_std['region']=df['region']
df_ensemble_std['name']=df['name']
df_ensemble_std['sector_type']=df['sector_type']
df_ensemble_std['ice_region_type']=df['ice_region_type']

df=df_ensemble_mean
dfstd=df_ensemble_std

for var in varx:
    df[var][((df['mmSLE_mean']<0)&(dfstd['mmSLE_mean'] >3))]=np.nan
dfstd.columns
#%% diagnostics

# plt.plot(df_ensemble_mean.mmSLE_2012)
# # plt.plot(df_ensemble_std.mmSLE_mean)
# # plt.plot(df_ensemble_mean.mmSLE_mean+df_ensemble_std.mmSLE_2012)
# # plt.plot(df_ensemble_mean.mmSLE_mean-df_ensemble_std.mmSLE_2012)

# # plt.scatter(df_ensemble_mean.mmSLE_mean,df_ensemble_mean.mmSLE_2012)

# #%%
# # plt.plot(df_ensemble_mean.valid_flag)

# plt.scatter(df.mmSLE_mean,dfstd.mmSLE_mean)

# df.name[df.mmSLE_mean<0]
# df.area[df.mmSLE_mean<0]
# df.mmSLE_mean[df.mmSLE_mean<0]
# df.name[dfstd.mmSLE_mean>3]
# df.mmSLE_mean[dfstd.mmSLE_mean>3]
# df.area[dfstd.mmSLE_mean>3]


#%% sector scenario and regional results

scenario_name_names=['TMB, all sectors','TMB, tidewater sectors','SMB, land terminating sectors','SMB, tidewater sectors']

scenario_names=['all_TMB','TW_TMB','LT_SMB']

# tot_area=np.sum(df1.area)
# tot_vol=np.sum(df1.volume)

if_icecaps=0
wo=0

for gamma_index in range(0,1):    
# for gamma_index in range(0,2):
    gamma=1.
    gamma_name='_gamma1.0'
    
    if gamma_index:
        gamma=1.25
        gamma_name='_gamma1.25'
    for scenario_name_index,scenario_name in enumerate(scenario_names):
        # if scenario_name_index>=0:
        if scenario_name_index==0: # all
        # if scenario_name_index==1: #TW
        # if scenario_name_index==2: #LT
            # for if_icecaps in range(0,1):
                if if_icecaps:
                    if_icecaps_name='_icecaps_only'
                    if_icecaps_name2=', ice caps only'
                    regions=['SW','CW','NW','NO','NE','CE'] ; n_regions=len(regions)
                else:
                    if_icecaps_name=''
                    if_icecaps_name2=''
                    regions=['SW','CW','NW','NO','NE','CE','SE','all'] ; n_regions=len(regions)
                regions=['all'] ; n_regions=len(regions)
                regions=['SW'] ; n_regions=len(regions)
                # for by_area in range(0,2):
                for by_area in range(1,2):
                    if by_area:
                        by_area_name='_by_area'
                    else:
                        by_area_name='_by_volume'
                    # cum=0.
                    volsum=0.
    
                if wo:
                    ofile='/Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/___'+\
                    scenario_name+\
                    '_AAR_vs_MB_all_sectors_regionally_ensemble'+\
                    if_icecaps_name+\
                    by_area_name+\
                    gamma_name+\
                    '.csv'
                
                    # ofile_all_years='/Users/jason/Dropbox/ELA/stats/imbalance_regional_sector_integration/___'+\
                    # scenario_name+\
                    # '_AAR_vs_MB_all_sectors_regionally_ensemble'+\
                    # if_icecaps_name+\
                    # by_area_name+\
                    # gamma_name+\
                    # '_all_years.csv'
                
                    out=open(ofile,'w+')
                    out.write(scenario_names[scenario_name_index]+'\n')
                    out.write('sector,area,area fraction,volume,volume fraction,imbalance perpetual 2000-2019,imbalance perpetual 2012,imbalance perpetual 2018,TMB\n')
                
                    ofile2='/Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/___'+scenario_name+\
                    '_AAR_vs_MB_all_sectors_regionally_ensemble'+\
                    if_icecaps_name+\
                    by_area_name+\
                    gamma_name+\
                    '_err.csv'
                
                    out_err=open(ofile2,'w+')
                    out_err.write('sector,area,area fraction,volume,volume fraction,imbalance perpetual 2000-2019,imbalance perpetual 2012,imbalance perpetual 2018,TMB\n')
                
                for region_index,region in enumerate(regions):
                    
                    # if region!='all':
                    #     v=df.region==region
                    # else:
                    #     v=df.region!='null'
                        
                    # tot_regional_area=np.sum(df.area[v])
                    # tot_regional_vol=np.sum(df.volume[v])

                    if scenario_name_index==0:
                        if region=='all':
                            tot_vol=sum(df.volume)
                            tot_area=sum(df.area)    
                            tot_regional_vol=sum(df.volume)
                            tot_regional_area=sum(df.area)
                            v=df.sector_type!='null'
                        else:
                            v=df.region==region
                            tot_vol=sum(df.volume)
                            tot_area=sum(df.area)
                            tot_regional_vol=sum(df.volume[v])
                            tot_regional_area=sum(df.area[v])                
                            print(region,tot_regional_area)
                            plt.plot(df.mmSLE_mean[v],'.')
                            print(sum(v),sum(df.mmSLE_mean[v]))

                    if scenario_name_index==1:
                        if region=='all':
                            tot_regional_vol=sum(df.volume[df.sector_type=='TW'])
                            tot_regional_area=sum(df.area[df.sector_type=='TW'])
                            v=df.sector_type=='TW'
                        else:
                            v=((df.region==region)&(df.sector_type=='TW'))
                            tot_regional_vol=np.nansum(df.volume[v])
                            tot_regional_area=sum(df.area[v])

                    if scenario_name_index==2:
                        if region=='all':
                            tot_regional_vol=sum(df.volume[df.sector_type=='LT'])
                            tot_regional_area=sum(df.area[df.sector_type=='LT'])
                            v=df.sector_type=='LT'
                        else:
                            v=((df.region==region)&(df.sector_type=='LT'))
                            tot_regional_vol=np.nansum(df.volume[v])
                            tot_regional_area=sum(df.area[v])

                    # if scenario_name_index==3:
                    #     if region=='all':
                    #         tot_regional_vol=sum(df.volume[df.sector_type=='TW'])
                    #         tot_regional_area=sum(df.area[df.sector_type=='TW'])
                    #     else:
                    #         v=((region==region)&(df.sector_type=='TW'))
                    #         tot_regional_vol=np.nansum(df.volume[v])
                    #         tot_regional_area=sum(df.area[v])

                    mm_SLE_mean=np.sum(df.mmSLE_mean[v])
                    mm_SLE_2012=np.sum(df.mmSLE_2012[v])
                    mm_SLE_2018=np.sum(df.mmSLE_2018[v])

                    mm_SLE_mean_err=np.sum(dfstd.mmSLE_mean[v])
                    mm_SLE_2012_err=np.sum(dfstd.mmSLE_2012[v])
                    mm_SLE_2018_err=np.sum(dfstd.mmSLE_2018[v])
                    
                    TMB=np.sum(df['mean mass flux'][v])
                    TMB_err=np.sum(dfstd['mean mass flux'][v])

                    print(region,mm_SLE_mean,tot_regional_area) #!!! integrate by area weighthing?
                    # print(np.sum(df.mmSLE_mean_truncated[v]))
                    # print(np.sum(df.mmSLE_mean[v]))
                    
                    # print(np.sum(df.mmSLE_mean[v])/np.sum(df.mmSLE_mean_truncated[v]))
                
                # #%% mean LT AAR
                # v=((df.sector_type=='LT')&(df.R>0.4)&(df.area>1000.)) # !!!
                # print('mean AAR0 LT',np.mean(df.AAR0[v]))
                
                    if wo:
                        out.write(region+\
                                      ','+str("%.0f"%tot_regional_area)+\
                                      ','+str("%.2f"%(tot_regional_area/tot_area))+\
                                      ','+str("%.0f"%tot_regional_vol)+\
                                      ','+str("%.2f"%(tot_regional_vol/tot_vol))+\
                                      ','+str("%.0f"%mm_SLE_mean)+\
                                      ','+str("%.0f"%mm_SLE_2012)+\
                                      ','+str("%.0f"%mm_SLE_2018)+\
                                      ','+str("%.0f"%TMB)+\
                                      ' \n')
                
                    # out_all_years.write(region+\
                    #               ','+str("%.0f"%tot_regional_area)+\
                    #               ','+str("%.2f"%(tot_regional_area/tot_area))+\
                    #               ','+str("%.0f"%tot_regional_vol)+\
                    #               ','+str("%.2f"%(tot_regional_vol/tot_vol))+\
                    #               ','+str("%.0f"%mm_SLE_mean)+\
                    #               ','+str("%.2f"%mm_SLE_x[0])+\
                    #               ','+str("%.2f"%mm_SLE_x[1])+\
                    #               ','+str("%.2f"%mm_SLE_x[2])+\
                    #               ','+str("%.2f"%mm_SLE_x[3])+\
                    #               ','+str("%.2f"%mm_SLE_x[4])+\
                    #               ','+str("%.2f"%mm_SLE_x[5])+\
                    #               ','+str("%.2f"%mm_SLE_x[6])+\
                    #               ','+str("%.2f"%mm_SLE_x[7])+\
                    #               ','+str("%.2f"%mm_SLE_x[8])+\
                    #               ','+str("%.2f"%mm_SLE_x[9])+\
                    #               ','+str("%.2f"%mm_SLE_x[10])+\
                    #               ','+str("%.2f"%mm_SLE_x[11])+\
                    #               ','+str("%.2f"%mm_SLE_x[12])+\
                    #               ','+str("%.2f"%mm_SLE_x[13])+\
                    #               ','+str("%.2f"%mm_SLE_x[14])+\
                    #               ','+str("%.2f"%mm_SLE_x[15])+\
                    #               ','+str("%.2f"%mm_SLE_x[16])+\
                    #               ','+str("%.2f"%mm_SLE_x[17])+\
                    #               ','+str("%.2f"%mm_SLE_x[18])+\
                    #               ','+str("%.2f"%mm_SLE_x[19])+\
                    #               ' \n')
                        # # pxxx
                    # print(region+\
                    #               ','+str("%.0f"%tot_regional_area)+\
                    #               ','+str("%.2f"%(tot_regional_area/tot_area))+\
                    #               ','+str("%.0f"%tot_regional_vol)+\
                    #               ','+str("%.2f"%(tot_regional_vol/tot_vol))+\
                    #               ','+str("%.0f"%mm_SLE_mean)+\
                    #               ','+str("%.0f"%mm_SLE_x[12])+\
                    #               ','+str("%.0f"%mm_SLE_x[18])+\
                    #               ','+str("%.0f"%spec_disequil_all_mean)+\
                    #               ','+str("%.0f"%spec_disequil_x[12])+\
                    #               ','+str("%.0f"%spec_disequil_x[18])+\
                    #               ','+str("%.0f"%(np.nansum(x)/20.))+\
                    #               ' \n')
                
                    if wo:
                        out_err.write(region+\
                                      ','+str("%.0f"%tot_regional_area)+\
                                      ','+str("%.2f"%(tot_regional_area/tot_area))+\
                                      ','+str("%.0f"%tot_regional_vol)+\
                                      ','+str("%.2f"%(tot_regional_vol/tot_vol))+\
                                      ','+str("%.0f"%mm_SLE_mean)+\
                                          '±'+str("%.0f"%(mm_SLE_mean_err))+\
                                      ','+str("%.0f"%mm_SLE_2012)+'±'+str("%.0f"%(mm_SLE_2012_err))+\
                                      ','+str("%.0f"%mm_SLE_2018)+'±'+str("%.0f"%(mm_SLE_2018_err))+\
                                      ','+str("%.0f"%TMB)+'±'+str("%.0f"%(TMB_err))+\
                                      ' \n')          
                
                if wo:
                    # out_all_years.close()
                    out.close()
                    out_err.close()
                msg='/Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/___*'+gamma_name+'*.csv'
                os.system('ls -lF '+msg)
                # os.system('cat /Users/jason/Dropbox/ELA/stats/imbalance_regional_sector_integration/*'+gamma_name+'*')
                os.system('cat '+msg+' > /Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/cat_AAR_vs_MB_together'+gamma_name+'.csv')