#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 10:29:25 2019

out to 

ELA/Figs/AAR_vs/
ELA/stats/

@author: jason

followed by /Users/jason/Dropbox/ELA/prog/integrate_results.py

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

ly='p'
wo=0
batch=0
prt=0
open_plot=0
plt_map=1
do_plot=1
do_effective_AAR0=1

if_icecaps_name=''
thresh=9

do_SMB=0
plt_name='SMB'
plt_name2='surface mass balance'

if do_SMB==0:
    plt_name='TMB'
    plt_name2='total mass balance'

def add_subplot_axes(ax,rect):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

# plt.rcParams['font.sans-serif'] = ['Georgia']
th=1 # line thickness
plt.rcParams["font.size"] = 12
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'

fn='/Users/jason/Dropbox/1km_grid2/sector_info_v3.csv'
#count,name,region,type,id,lat,lon,area
#os.system('ls -lF '+fn)
#os.system('head '+fn) 
df = pd.read_csv(fn, delimiter=",")
# df.columns
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
RCM_index0=0 ; RCM_index1=1 # MAR

if batch:
    SID_index0=0 ; SID_index1=3
    ALB_index0=0 ; ALB_index1=3
    volume_index0=0 ; volume_index1=1
    RCM_index0=0 ; RCM_index1=2

# ------------------------------------------------ loop over SID uncertainty
devname=["lo","mid","hi","const"]

SID_devname=["0.9","1.0","1.1","const"]

ccc=0
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

                if wo:
                    out_fn='/Users/jason/Dropbox/ELA/stats/imbalance/sector_scale/'+\
                        'ALB'+devname[ALB_dev_index]+'_'+RCM[mm]+'_'+\
                            volume_name[volume_index]+'_SID'+SID_devname[SID_dev_index]+\
                                if_icecaps_name
                    out=open(out_fn+'.csv','w')
                    out.write('name,region,sector_type,ice_region_type,lat,lon,volume,area,mean mass flux,AAR_mean_truncated,AAR0_truncated,alpha_mean_truncated,mmSLE_mean_truncated,AAR_mean,AAR0,alpha_mean,mmSLE_mean,AAR_2012,alpha_2012,mmSLE_2012,AAR_2018,alpha_2018,mmSLE_2018,R,N,SMB_mean,SID_mean,mm SLE LT_effective ÷ mm SLE ,valid_flag\n')

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + loop over sectors
                AARs0=np.zeros((n_large_enough_sectors,n_years))*np.nan
                AARs1=np.zeros((n_large_enough_sectors,n_years))*np.nan
                MBs=np.zeros((n_large_enough_sectors,n_years))*np.nan
                
                ccc+=1
                print(ccc)
                for k in range(n_large_enough_sectors):
                    # if df.name[k]!='null':
                    # if k>470:
                    # if k==4
                    # if name[k]=='OSTENFELD_GLETSCHER':
                    # # if region[k]=='NO':
                    # if df.name[k]=='STORSTROMMEN': 
                    # if df.name[k]=='UKAASORSUAQ': 
                    # if df.name[k]=='AB_DRACHMANN_GLETSCHER_L_BISTRUP_BRAE': 
                    if df.name[k]=='JAKOBSHAVN_ISBRAE':
                    # if df.name[k]=='SAQQAP-MAJORQAQ-SOUTHTERRUSSEL_SOUTHQUARUSSEL':

                        # print(k,df.name[k])
                        if prt:
                            print()
                            print(k,df.name[k])
    
                        # print(name[k])
                        # =================================== ELA
                        fn='/Users/jason/Dropbox/ELA/stats/TSL/'+namex[k]+'_ELA_v'+versionx+'_'+devname[ALB_dev_index]+'.csv'
                        df_ELA = pd.read_csv(fn)

                        # =================================== SMB
                        fn='/Users/jason/Dropbox/RCM/stats/'+df.name[k]+'_v'+versionx+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SMB.csv'
                        # https://www.dropbox.com/sh/i6uwag0uod54ngg/AADrDk24ndA-tpwOnPZOeUzsa?dl=0
                        df_SMB = pd.read_csv(fn, delimiter=",")
                        if df["sector_type"][k]=='TW':
                            # =================================== SID
                            fn='/Users/jason/Dropbox/SID/Mankoff_2019/polygons/'+\
                            name[k]+'_SID_ALB'+devname[ALB_dev_index]+'_SID'+SID_devname[SID_dev_index]+\
                                '_v'+versionx+'.csv'
                            # https://www.dropbox.com/sh/j08ic8tmffm4rpp/AAA0YGxdh9XYoK91NXQfHHR6a?dl=0
                            df_SID = pd.read_csv(fn)
                                                        
                        if df["sector_type"][k]=='TW': 
                            MB=df_SMB.smb-df_SID.SID
                        if df["sector_type"][k]=='LT':
                            MB=df_SMB.smb
                        if df["sector_type"][k]=='TW':
                            SID_mean=np.nanmean(df_SID.SID)
                        else:
                            SID_mean=0.

                        MB_mean=np.nanmean(MB)
                        SMB_mean=np.nanmean(df_SMB.smb)

                        # x=MB
                        # y=df_ELA.AAR
    
                        # idx=sectorx.values[k]
                        # sector=str('%03d'%idx)
                                 
                        year = df_ELA["year"]
                        AAR = df_ELA.AAR
                        # ELA = df_ELA["ELA"]
                        
                        if df.name[k]=='STORSTROMMEN': AAR[1]=0.
                        # df.loc[df[((df.name[k]=='STORSTROMMEN')&()), 'B'] = new_val

                        if df.name[k]=='AB_DRACHMANN_GLETSCHER_L_BISTRUP_BRAE': AAR[1]=0.
                        if df.name[k]=='HAGEN_BRAE': AAR[16]=0.
                                            
                        v=np.where((df_ELA.AAR > 0.)&(df_ELA.AAR < 1.)) ; y=df_ELA.AAR[v[0]].values ; x=MB[v[0]].values ; year=year[v[0]].values
                        v=v[0]

                        # populate arrays to output
                        AARs0[k,v]=y
                        if df["sector_type"][k]=='LT':AARs1[k,v]=y
                        MBs[k,:]=MB.values
                        
                        c=len(v)
                        if c<thresh:
                            lost_count+=1
                        alpha_mean=0.
                        AAR_2012=0.

                        if c>=thresh:
                            statsx=stats.pearsonr(x, y)
                            # plt.scatter(x,y)
                            R=statsx[0]
                            Rthresh=0.
                            if R<Rthresh:
                                lost_count+=1
                                # print("invalid",name[k],region[k],area[k]) # !!!!!!!!!!
                            if R>Rthresh:
                #                print(k,c,name[k],region[k],gl_type[k],'R',R)
                                valid_flag[k]=1
                                valid_vol_sum+=volume[k]
                                valid_area_sum+=area[k]
                            
                                if do_plot==1:
                                    plt.close()
                                    fig, ax = plt.subplots(figsize=(6,4))

                                    if do_effective_AAR0==0:ax.scatter(x,y,color='w',s=1)#,size=1)
                                    ax.set_xlabel(plt_name2+', Gt/y')
                                    ax.set_ylabel('AAR')
                                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

                                coef = np.polyfit(x,y,1)
                                m,AAR0_truncated = np.polyfit(x, y, 1)
                                # poly1d_fn = np.poly1d(coef) 
                        
                                xx = np.linspace(np.min(x), np.max(x),2)
                                yx = m*xx+AAR0_truncated
                                
                                if do_plot==1:
                                    ax.plot(xx, yx, '--',color='grey',linewidth=th*0.7)
                                
                                AAR_mean_truncated=np.mean(y)
                                mean_mass_flux=np.mean(x)
                                alpha_mean_truncated=AAR_mean_truncated/AAR0_truncated
                #                if alpha_mean<0.5:
                ##                    print('warning',name2[k],str("%9.5f"%AAR_mean).lstrip(),str("%9.5f"%AAR0).lstrip(),str("%9.5f"%alpha_mean).lstrip())
                #                    print('warning',name2[k],str("%9.5f"%alpha_mean).lstrip())
                                
                                # if area[k]<900:alpha_mean=0.8
                                # imbalance_mean=(1-alpha_mean)*volume[k]
                                # mm_SLE_mean=imbalance_mean/362.
                                gamma=1
                                imbalance_mean_truncated=(alpha_mean_truncated**gamma-1)*df["volume"][k]
                                mm_SLE_mean_truncated=-imbalance_mean_truncated/362.
                                
                                AAR_2012=0.
                                AAR_2018=0.
                                for kk in range(0,c):
                                    if year[kk]==2012:AAR_2012=AAR[kk]
                                    if year[kk]==2018:AAR_2018=AAR[kk]
                                    
                                alpha_2012=AAR_2012/AAR0_truncated
                                imbalance_2012=(alpha_2012**gamma-1)*df["volume"][k]
                                mm_SLE_2012=-imbalance_2012/362.

                                alpha_2018=AAR_2018/AAR0_truncated
                                imbalance_2018=(alpha_2018**gamma-1)*df["volume"][k]
                                mm_SLE_2018=-imbalance_2018/362.

#                                 print(AAR0_truncated)
# #%%
# ----------------------------------------------------------------------------- effective ablation area
                                if ( (df["sector_type"][k]=='TW') & (do_effective_AAR0) ):
                                    AAR0_target=0.782 #!!!
                                    AAR_anoms=df_ELA.AAR[v].values-AAR_mean_truncated
                                    std_anoms=AAR_anoms/np.std(AAR_anoms)
                                    # plt.close()
                                    # plt.plot(AAR_anoms/np.std(AAR_anoms),'o')
                                    
                                    # std_AAR_anoms=np.std(AAR_anoms)/AAR_mean
                                    
                                    # AAR_adjusted=df_ELA.AAR/(AAR_mean/AAR0_target)
                                    # y=AAR_adjusted
                                    area_accum=AAR_mean_truncated*df["area"][k]
                                    area_abl=df["area"][k]*(1-AAR_mean_truncated)
                                    # asas

                                    area_fraction=0.001*df["area"][k]
                                    vals=[]
                                    for i in range(0,800):
                                        area_abl_adj=area_abl+i*area_fraction
                                        area_adj=area_abl_adj+area_accum
                                        AAR_adj=area_accum/area_adj
                                        y=AAR_anoms+AAR_adj
                                        m,AAR0 = np.polyfit(x, y, 1)
                                        # vals.append(abs(AAR0_truncated-AAR0_target))
                                        if abs(AAR0-AAR0_target)<0.01:break
                                        # print(i,AAR_adj,AAR0,abs(AAR0-AAR0_target))

                                    if do_plot:
                                       # ax.scatter(x, y, facecolors=None, edgecolors=None,zorder=9)
                                       for i in range(0,20):
                        #                    print(i,SMB[i], AAR[i])
                                        if ( (y[i] > 0.0) & (y[i] < 1) ):
                                            ax.text(x[i], y[i], str(year[i])[2:4], \
                                                 horizontalalignment='center',verticalalignment='center', \
                                                 color='r')
                                    
                                    statsx=stats.pearsonr(x, y)
                                    R=statsx[0]
                                    oneminusp=1-statsx[1]
                                    oneminusp_string=' = '+str("%8.3f"%oneminusp).lstrip()
                                    
                                    # populate arrays to output
                                    AARs1[k,v]=y
                                    annotate=1
                                    if round(oneminusp, 3)==1.0:
                                        oneminusp_string='>0.999'
                                    # print(round(oneminusp, 3))
                                    
                                    coef = np.polyfit(x,y,1)
                                    m,AAR0 = np.polyfit(x, y, 1)
                                    # poly1d_fn = np.poly1d(coef) 
                                    
                                    xx = np.linspace(np.min(x), 0,2)
                                    yx = m*xx+AAR0
                                    
                                    detrended=y-(m*x+AAR0)
                                    std_resid=np.nanstd(detrended)
                                    
                                    if do_plot:
                                        plt.title(df.name[k])
                                        plt.plot(xx, yx, '--',color='grey',linewidth=th*2)
                                                             
                                    gamma=1
                                    
                                    AAR_mean=np.mean(y)
                                    alpha_mean=AAR_mean/AAR0
                                    dA=-(alpha_mean-1)*df["area"][k] #!!! OK?
                                    imbalance_mean=(alpha_mean**gamma-1)*df["volume"][k]
                                    mm_SLE_mean=-imbalance_mean/362.
                                    
                                    for kk in range(0,c):
                                        if year[kk]==2012:AAR_2012=AAR[kk]
                                        if year[kk]==2018:AAR_2018=AAR[kk]
                                        if prt:print(year[kk],x[kk],y[kk])
                                    alpha_2012=AAR_2012/AAR0_truncated
                                    imbalance_2012=(alpha_2012**gamma-1)*df["volume"][k]
                                    mm_SLE_2012=-imbalance_2012/362.
    
                                    alpha_2018=AAR_2018/AAR0_truncated
                                    imbalance_2018=(alpha_2018**gamma-1)*df["volume"][k]
                                    mm_SLE_2018=-imbalance_2018/362.
                                    if prt:
                                        print("committed area loss = {:.0f}".format(dA)+"±{:.0f}".format((dA*0.28)))
                                        print("imbalance = {:.0f}".format(imbalance_mean)+' Gt')
                                        print("mm SLR = {:.1f}".format(mm_SLE_mean))
                                else:
                                    AAR_mean=AAR_mean_truncated
                                    alpha_mean=alpha_mean_truncated
                                    mm_SLE_mean=mm_SLE_mean_truncated
                                    AAR0=AAR0_truncated
                                    
# ----------------------------------------------------------------------------- end effective ablation area

                                if do_plot==1:
                                    xos=abs(np.nanmean(x))*0.06
                                    xx = np.linspace(np.min(x)-xos, 0,2)
                                    yx = [AAR0,AAR0]
                                    ax.plot(xx, yx, '--', color='b',linewidth=th*0.7) # blue line
                #                    ax.text(np.min(x)+xos,yx[0]+yx[0]*0.005,'AAR$_0$ = '+str("%8.3f"%AAR0).lstrip(),color='b')
                                    ax.text(np.min(x)+xos,AAR0,'AAR$_0$ = '+str("%.3f"%AAR0).lstrip(),color='b')
                        
                                    ymin=0.5
                                    yos=0.001
                                    if np.min(y)-yos < ymin:ymin=0.
                    
                                    xx = [0,0]
                    #                yx = [np.min(y),AAR0]
                                    yx = [ymin,1.]
                                    ax.plot(xx, yx, '--', color='orange',linewidth=th*0.7)
                                
                                    ymax=np.max(y)+yos
                                    if AAR0>ymax:ymax=AAR0+yos
                                    ax.set_xlim(np.min(x)-xos, np.max(x)+xos)
                                    ax.set_ylim(np.min(y)-yos, ymax)
                                                
                                    ax.set_title(str(name[k].title()+'\n'+gl_type[k].replace('_',' ')+\
                                                     ' '+df["sector_type"][k]+' '+'sector, '+"%4.1f"%lat[k])+' N, '+\
                                                 str("%5.1f"%lon[k])+', '+f'{area[k]:n}'+' $km^2$')
                                    
                                    # 2 digit year
                                    if do_effective_AAR0==0:
                                        for i in range(0,20):
                        #                    print(i,SMB[i], AAR[i])
                                            if ( (AAR[i] > 0.0) & (AAR[i] < 1) ):
                                                ax.text(MB[i], AAR[i], str(df_ELA["year"][i])[2:4], \
                                                     horizontalalignment='center',verticalalignment='center', \
                                                     color='r')
                                                
                                    


                                    #-----------------------------------------------------------------
                                    # cc=0. ; dy=44 ; xx0=1660 ; yy0=490 ; fs=9
                                    cc=0. ; dy=44 ; xx0=1450 ; yy0=490 ; fs=9
                                    #-----------------------------------------------------------------
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on map area
                                        xx0, yy0-cc*dy, 'R = '+str("%.3f"%R),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, '$\overline{AAR}$ = '+str("%.3f"%(AAR_mean)),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    aa=r'$\alpha$'
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, aa+' = '+str("%.4f"%(alpha_mean)),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, '$\overline{SMB}$: '+str("%.2f"%(SMB_mean)+' Gt/y'),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, '$\overline{SID}$: '+str("%.2f"%(SID_mean)+' Gt/y'),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, '$\overline{MB}$: '+str("%.2f"%(MB_mean)+' Gt/y'),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, 'volume: '+str("%.2f"%(100*volume[k]/tot_vol)+' %'),
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, 'N = '+str(c)+' years',
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.
                                    #-----------------------------------------------------------------
                                    ax.text(  # position text absolutely at specific pixel on image
                                        xx0, yy0-cc*dy, 'imbalance:\n'+str("%8.3f"%mm_SLE_mean).lstrip()+' mm SLE',
                                        ha='left', va='top',
                                        transform=None,color='g', fontsize=fs)
                                    cc+=1.                
                                    cc+=1.                
                                    #-----------------------------------------------------------------
                                    # ax.text(  # position text absolutely at specific pixel on image
                                    #     xx0, yy0-cc*dy, 'imbalance 2012:\n'+str("%8.3f"%mm_SLE_2012).lstrip()+' mm SLE',
                                    #     ha='left', va='top',
                                    #     transform=None,color='g', fontsize=fs)
                                    # cc+=1.                
                                    # cc+=1.
                                    #-----------------------------------------------------------------
                                    # ax.text(  # position text absolutely at specific pixel on image
                                    #     xx0, yy0-cc*dy, 'alb selection:'+devname[ALB_dev_index],
                                    #     ha='left', va='top',
                                    #     transform=None,color='g', fontsize=fs)
             
                                    #-----------------------------------------------------------------
                #                    ax.text(  # position text absolutely at specific pixel on image
                #                        xx0, yy0-cc*dy, '2 digit year',
                #                        ha='left', va='top',
                #                        transform=None,color='r', fontsize=fs)
                #                    cc+=1.
                                    props = dict(boxstyle='round', facecolor='w',edgecolor='grey',alpha=0.6)
                    #                    ax.text(0.8, 0.15, 'J.Box, GEUS', transform=ax.transAxes, fontsize=9,
                    #                            verticalalignment='top', bbox=props)
                    #                    ax.text(0.7, 0.17,'AAR$_0$ = '+str("%8.4f"%AAR0).lstrip(),color='b',
                    #                            transform=ax.transAxes,#, fontsize=9
                    #                            verticalalignment='top', bbox=props)
                                    ax.text(0.77, 0.09,'2-digit year',color='r',
                                            transform=ax.transAxes,#, fontsize=9
                                            verticalalignment='top', bbox=props)

                                    if plt_map:
                                        ax = fig.add_subplot(111)
                                        rect = [0.79,0.48,0.7,0.7]
                                        
                                        ax1 = add_subplot_axes(ax,rect)
                                        fn='/Users/jason/Dropbox/ELA/ancil/mouginot/Mouginot_2019/Greenland_Basins_PS_v1.4.2.shp'
                                        GL = gpd.read_file(fn)
                                        temp=GL.NAME
                                        for l in range(len(GL)):
                                            if l!=df["sector_id"][k]:
                                                temp[l]='t'
                                        GL.plot(temp,ax=ax1,cmap='viridis_r')

                                        ax1.axis('off')
                                        ax.axis('off')
                
                        
                                if do_plot==1:
                                    if ly == 'x':
                                        plt.show()
                                    
                                    if ly == 'p':
                                        figpath='/Users/jason/Dropbox/ELA/Figs/AAR_vs_'+plt_name+'/'
                                        # os.system('mkdir -p '+figpath)
                                        # figpath_regional='/Users/jason/Dropbox/ELA/Figs/AAR_vs_'+plt_name+'/_by_region/'+region[k]+'/'
                                        # os.system('mkdir -p '+figpath_regional)
                    
                                        #----------------------------------------------- by name
                                        figpath_regional=figpath+'_by_name/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)                                        # plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        # #----------------------------------------------- by lat                                            
                                        # figname=figpath_regional+str("%6.2f"%lat[k]).lstrip()+'_'+region[k]+'_'+name[k]+'_AAR_vs_'+plt_name+'_'+versionx+'_'+devname[dev_index]+'.png'
                                        # plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        # figpathx=figpath+'_by_latitude/' ; os.system('mkdir -p '+figpathx)
                                        # os.system('/bin/cp '+figname+' '+figpathx)
                                        #----------------------------------------------- AAR0
                                        figpath_regional=figpath+'_by_AAR0/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+str("%.3f"%AAR0)+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        #----------------------------------------------- R
                                        figpath_regional=figpath+'_by_R/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+str("%.3f"%R)+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        #----------------------------------------------- Region
                                        figpath_regional=figpath+'_by_region/' ; os.system('mkdir -p '+figpath_regional)
                                        figpath_regional=figpath+'_by_region/'+region[k]+'/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+str("%.2f"%(100*volume[k]/tot_vol))+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)

                                        # #----------------------------------------------- by alpha
                                        # figname=figpath_regional+str("%8.5f"%(alpha_mean)).lstrip()+'_'+region[k]+'_'+name[k]+'_AAR_vs_'+plt_name+'_'+versionx+'_'+devname[dev_index]+'.png'
                                        # plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        # figpathx=figpath+'_by_alpha/' ; os.system('mkdir -p '+figpathx)
                                        # os.system('/bin/cp '+figname+' '+figpathx)
                                        # #----------------------------------------------- by type
                                        # figname=figpath_regional+gl_type[k]+'_'+str("%8.3f"%R).lstrip()+'_'+region[k]+'_'+name[k]+'_AAR_vs_'+plt_name+'_'+versionx+'_'+devname[dev_index]+'.png'
                                        # plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        # figpathx=figpath+'_by_type/' ; os.system('mkdir -p '+figpathx)
                                        # os.system('/bin/cp '+figname+' '+figpathx)
                                        #----------------------------------------------- by volume fraction
                                        figpath_regional=figpath+'_by_volfrac/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+str("%.2f"%(100*volume[k]/tot_vol))+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        # #----------------------------------------------- by area
                                        # figpath_regional=figpath+'_by_area/' ; os.system('mkdir -p '+figpath_regional)
                                        # figname=figpath_regional+str("%.2f"%(100*area[k])+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        # plt.savefig(figname, bbox_inches='tight', dpi=250)

                                        #----------------------------------------------- by SLR
                                        figpath_regional=figpath+'_by_vSLR/' ; os.system('mkdir -p '+figpath_regional)
                                        figname=figpath_regional+str("%.2f"%(mm_SLE_mean))+'_'+name[k]+'_'+devname[ALB_dev_index]+'_'+RCM[mm]+'_SID'+SID_devname[SID_dev_index]+'.png'
                                        plt.savefig(figname, bbox_inches='tight', dpi=250)
                                        
                                        if open_plot:
                                            os.system('open '+figname)
                        #                os.system('ls -lF '+figname)
                        
                        if wo:
                            out.write(name[k]+','+\
                                          region[k]+','+\
                                          df["sector_type"][k]+','+\
                                          gl_type[k]+','+\
                                          str("%.2f"%lat[k])+','+\
                                          str("%.2f"%lon[k])+','+\
                                          str("%.0f"%volume[k])+','+\
                                          str("%.1f"%area[k])+','+\
                                          str("%.1f"%mean_mass_flux)+','+\
                                          str("%.5f"%AAR_mean_truncated)+','+\
                                          str("%.5f"%AAR0_truncated)+','+\
                                          str("%.5f"%alpha_mean_truncated)+','+\
                                          str("%.4f"%mm_SLE_mean_truncated)+','+\
                                          str("%.5f"%AAR_mean)+','+\
                                          str("%.5f"%AAR0)+','+\
                                          str("%.5f"%alpha_mean)+','+\
                                          str("%.4f"%mm_SLE_mean)+','+\
                                          str("%.5f"%AAR_2012)+','+\
                                          str("%.5f"%alpha_2012)+','+\
                                          str("%.4f"%mm_SLE_2012)+','+\
                                          str("%.5f"%AAR_2018)+','+\
                                          str("%.5f"%alpha_2018)+','+\
                                          str("%.4f"%mm_SLE_2018)+','+\
                                              str("%.3f"%R).lstrip()+','+\
                                          str(c)+','+\
                                          str("%.2f"%SMB_mean)+','+\
                                          str("%.2f"%SID_mean)+','+\
                                          str("%.2f"%(mm_SLE_mean/mm_SLE_mean_truncated))+','+\
                                          # str("%6.1f"%imbalance_2012)+','+\
                                          # str("%6.1f"%mm_SLE_2012)+','+\
                                          # sector+','+\
                                          str("%.0f"%valid_flag[k])+
                                          ' \n')
           
                if wo:
                    out.close()
                    # AARs.astype('float32').tofile(out_fn+'2x473x20.AARs.npy')
                    # MBs
                    # MBs.astype('float32').tofile(out_fn+'1x473x20.MBs.npy')
                    pd.DataFrame(MBs).to_csv(out_fn+'.MBs.csv')
                    pd.DataFrame(AARs0).to_csv(out_fn+'.AARs0.csv')
                    pd.DataFrame(AARs1).to_csv(out_fn+'.AARs1.csv')
                #    os.system('cat '+out_fn)
                    # os.system('wc -l '+out_fn+'.csv')
                    # os.system('cat '+out_fn+'.csv')
                            
                        # print("N with insufficient sample",lost_count)
                        # print("N with sufficient sample",n_sectors-lost_count)
                        # print("valid_vol frac",valid_vol_sum/tot_vol)
                        # print("valid_area frac",valid_area_sum/tot_areal)
                        
                        # valid_flag=np.asarray(valid_flag)
                        
                        # v=np.where(valid_flag<1)
                        # print("area invalid",sum(area[v[0]]))
                        # print("mean area invalid",np.mean(area[v[0]]))
                        
                        # v=np.where(valid_flag>0)
                        # print("area valid",sum(area[v[0]]))
                        # print("mean area valid",np.mean(area[v[0]]))
                        # print("imbalance_mean",imbalance_mean)
                        # # print("imbalance_2012",imbalance_2012)
                        # print("alpha_mean",alpha_mean)
                        # print("alpha_2012",alpha_2012)