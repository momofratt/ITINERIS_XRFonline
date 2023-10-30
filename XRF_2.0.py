#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:19:55 2023

@author: cosimo
"""

from read_gupix_output_2 import get_gupix_tables
import matplotlib.pyplot as plt
import pandas as pd
from XRF_functions_2 import plot_concentrations, get_std_2fwhm, merge_and_plot
from numpy import linspace, nan

markers = ['o','v','*','s','+','x']
colors = ['#de1d90','#a52a2a','#4682b4','#6b8e23','#f8766d']
colors = ['black','grey','lightgrey']


##################################################################################################
###################################      PARAMETRI MISURA     ####################################
# STRAS Terni He
# filenm_prefix      = 'terni_20231024_He_STD0'
# pixe_filenm        = 'conc_STRAS_Terni.csv'
# pixe_err_filenm    = 'conc_STRAS_Terni_errors.csv' # file contenente errori PIXE
# filter_prefix      = 'stras'
# res_folder         = './results/20231027/He/'
# pixe_nan           = 'BL' # valore nan nel file di risultati pixe
# nfiles             = 14
# time_str_filter    = '1800s' # time string reported in the filename of measured filters 
# voltage_str_legend = '38kV'  # stringa usata per legenda
# voltage_str_infile = ''      # stringa usata nei nomi file delle misure di standard e teflon
# min_std            = 5       # minuti di misura standard
# min_sample         = 30      # minuti di misura campioni
# I                  = 0.1     # corrente in mA
# T                  = min_std*60  # tempo di misura degli standard in secondi
# err_perc_standard  = 0.05        # errore percentuale sugli standard
# mylar_prefix       = 'standard_blk_'
# filter_blk_prefix  = 'stras_blk_'
# filtri             = ['stras22','stras45']

# TEFLON WHAT He
filenm_prefix      = 'terni_20231024_He_STD0'
pixe_filenm        = 'conc_PIXE_teflon_watman.csv'
pixe_err_filenm    = 'conc_PIXE_teflon_watman_errors.csv' # file contenente errori PIXE
filter_prefix      = 'T'
res_folder         = './results/20231027/He/'
pixe_nan           = 'BL' # valore nan nel file di risultati pixe
nfiles             = 14
time_str_filter    = '1800s' # time string reported in the filename of measured filters 
voltage_str_legend = '38kV'  # stringa usata per legenda
voltage_str_infile = ''      # stringa usata nei nomi file delle misure di standard e teflon
min_std            = 5       # minuti di misura standard
min_sample         = 30      # minuti di misura campioni
I                  = 0.1     # corrente in mA
T                  = min_std*60  # tempo di misura degli standard in secondi
err_perc_standard  = 0.05        # errore percentuale sugli standard
mylar_prefix       = 'standard_blk_'
filter_blk_prefix  = 'teflon_blk_'
filtri             = ['T3','T5']

# STRAS Terni
# filenm_prefix      = 'terni_20231025_STD0'
# pixe_filenm        = 'conc_STRAS_Terni.csv'
# pixe_err_filenm    = 'conc_STRAS_Terni_errors.csv' # file contenente errori PIXE
# filter_prefix='stras'
# res_folder = './results/20231027/'
# pixe_nan ='BL' # valore nan nel file di risultati pixe
# nfiles             = 25
# time_str_filter    = '1800s' # time string reported in the filename of measured filters 
# voltage_str_legend = '38kV'  # stringa usata per legenda
# voltage_str_infile = ''      # stringa usata nei nomi file delle misure di standard e teflon
# min_std            = 5       # minuti di misura standard
# min_sample         = 30      # minuti di misura campioni
# I                  = 0.1     # corrente in mA
# T                  = min_std*60  # tempo di misura degli standard in secondi
# err_perc_standard  = 0.05        # errore percentuale sugli standard
# mylar_prefix       = 'standard_blk_'
# filter_blk_prefix  = 'stras_blk_'
# filtri             = ['stras22','stras24','stras44','stras66','stras67']

# teflon whatman
# filenm_prefix = 'tef_what_20231024_STD0'
# pixe_filenm = 'conc_PIXE_teflon_watman.csv' # file contenente concentrazioni PIXE
# pixe_err_filenm = 'conc_PIXE_teflon_watman_errors.csv' # file contenente errori PIXE
# filter_prefix='T' # prefisso nome filtri

# nfiles = 20
# time_str_filter = '1800s' # time string reported in the filename of measured filters 
# voltage_str_legend = '38kV' # stringa usata per legenda
# voltage_str_infile = '' # stringa usata nei nomi file delle misure di standard e teflon
# min_std = 5 # minuti di misura standard
# min_sample = 30 # minuti di misura campioni
# I=0.1 # corrente in mA
# T=min_std*60 # tempo di misura degli standard in secondi
# err_perc_standard=0.05 # errore percentuale sugli standard
# mylar_prefix = 'standard_blk_'
# filter_blk_prefix = 'teflon_blk_'
# filtri = ['T3','T5']

# misura luglio
# filenm_prefix = 'teflon_what_STD0'
# nfiles = 28
# time_str_filter = '' # time string reported in the filename of measured filters 
# voltage_str_legend = '39kV' # stringa usata per legenda
# voltage_str_infile = '39kV' # stringa usata nei nomi file delle misure di standard e teflon
# min_sample = 5 # minuti di misura campioni
# min_std = 2 # minuti di misura standard
# I=0.1 # corrente in mA
# T=min_std*60 # tempo di misura degli standard in secondi
# err_perc_standard=0.05 # errore percentuale sugli standard
# mylar_prefix = 'std_blk_'
# filtri = ['T3', 'T4', 'T5','T9']
###################################################################################################
#####################################  TEFLON-WHATMAN   ###########################################
# cerca files generati da batch gupix e rinominali con il nome del file dello spettro
files_list=[]
for i in range(1,nfiles+1):
    try:
        file = open('results/' + filenm_prefix + str(i).zfill(2)+'.txt')
    except:
        continue
    
    line = file.readline()
    while 'File:' not in line:
        line = file.readline()
    filenm = line.split()[0].replace('File:','').split('_re')[0].split('.')[0].replace('_gupix','')
    files_list.append([filenm_prefix + str(i).zfill(2)+'.txt',filenm])

files_list.sort(key = lambda row: row[1]) # ordina secondo i nome file

for i in range(len(files_list)):  # rinomina standard con prefisso std_
    if not ((files_list[i][1].startswith(filter_prefix)) & (files_list[i][1][1] != 'i' ) & (files_list[i][1][1] != 'e' )):
        files_list[i][1] = 'std_' + files_list[i][1]
    
#######################################
########### # plot 2-FWHM  ############
#######################################

for filter_num in filtri:
    filters=[]
    fig, ax = plt.subplots(1,1, figsize=(8,10), sharey=True, sharex=True)
    colors     = ['blue','violet','green','red','orange','gray','salmon','turquoise','slateblue','seagreen']
    
    for i in range(len(files_list)):
        if (files_list[i][1].startswith(filter_num)) &  (files_list[i][1][1] != 'i' ) & (files_list[i][1][1] != 'e' ) & (time_str_filter in files_list[i][1]):
            filters.append( files_list[i] )
    
        
    for file, i in zip(filters, range(len(filters))):
        
        t1, t2, t3 = get_gupix_tables('./results/'+file[0])
    
        ax, t2 = merge_and_plot(t2, t3, ax, file[1], True, colors[i], '_', 
                                text=False, LOD=True, LOD_color=True, custom_lab=file[1])
        # fig.tight_layout()
        
    fig.savefig(res_folder+'2fwhm+lod_'  + filenm_prefix.replace('_STD0','') + '_' + filter_num+'_1-39kv.png', format='png', dpi=300)


##############   CALIBRAZIONE

standard_39kv   = [f for f in files_list if ('std_' in f[1]) 
                     & ('blk' not in f[1]) & (voltage_str_infile in f[1])]
mylar_39kv      = [f for f in files_list if (mylar_prefix in f[1]) ]
teflon_blk_39kv = [f for f in files_list if (filter_blk_prefix in f[1]) ]

# std_conc_39kv =  pd.DataFrame(data={'ele':[ 'Al', 'Na',  'K', 'Ca', 'Si', 'Mg', 'Cu', 'Fe', 'Pb', 'Zn', 'Mn', 'Ti',  'S', 'Cr', 'Ni', 'Mo'   ],
#                                    'conc':[ 40.8, 18.3, 26.0, 28.7, 38.4, 21.6, 16.7, 24.0, 43.4, 15.4, 24.2, 26.2, 11.3,   18.9  , 13.6   , 32.0   ],
#                                    'line':[ 'K',  'K',   'K',  'K',  'K',  'K',  'K',  'K',  'LA',  'K',  'K',  'K',  'K', 'K', 'K' , 'K' ]}) # valori di concentrazione elementale sui filtri misurati

std_conc_39kv =  pd.DataFrame(data={'ele':[ 'Al', 'Na',  'K', 'Ca', 'Si', 'Mg',   'S' ],
                                   'conc':[ 40.8, 18.3, 26.0, 28.7, 38.4, 21.6, 11.3   ],
                                   'line':[ 'K',  'K',   'K',  'K',  'K',  'K',  'K' ]}) # valori di concentrazione elementale sui filtri misurati


##########################################
############## sensitivity ###############
##########################################

peak_area_std_39kv = get_std_2fwhm(standard_39kv, mylar_39kv, std_conc_39kv, area = 'Peak_Area')


peak_area_std_39kv.insert(len(peak_area_std_39kv.columns), 'sensitivity', 0)


for ele in peak_area_std_39kv['element']:
    peak_area_std_39kv.loc[peak_area_std_39kv['element']==ele, 'sensitivity'] = \
                                (peak_area_std_39kv['peak_area'] - peak_area_std_39kv['mylar_peak_area'])/I/T/ \
                                    std_conc_39kv[std_conc_39kv['ele']==ele]['conc'].iat[0]

peak_area_std_39kv = peak_area_std_39kv.sort_values('Z')


peak_area_all = peak_area_std_39kv
peak_area_all = peak_area_all.sort_values('Z')

fig, ax = plt.subplots(1,1, figsize =(6,4))

ax.scatter(peak_area_all['element'], peak_area_all['sensitivity'], c='black', marker='*',label = voltage_str_legend,zorder=2)
ax.grid(which='both')
ax.legend()
ax.set_ylabel('sensitivity [cts/mA/($\mu$g/cm$^2$)]')
ax.set_yscale('log')
# ax.set_ylim([1,1000])
plt.savefig(res_folder+ filenm_prefix.replace('_STD0','') + '_sensitivity_Rh_anode.png', format='png', dpi=300)
                   
print("\n\nmylar blank contribution to sensitivity [%]:\n",round(1/(peak_area_std_39kv.set_index('element')['peak_area']/ peak_area_std_39kv.set_index('element')['mylar_peak_area'])*100,2))

#################################################
############### confronto XRF-PIXE ##############
################################################# 

mdl_frame_39kv = pd.DataFrame()

for i, fil in zip(range(len(filtri)),filtri):
    
    fig, ax= plt.subplots(1,1, figsize=(8,6))

    ############### plot XRF ###################

    filters_39kv=[]

    for i in range(len(files_list)):
        if (files_list[i][1].startswith(fil)) & \
                    ((files_list[i][1][1] != 'i' ) | (files_list[i][1][1] != 'e' )) & \
                    (time_str_filter in files_list[i][1]) & \
                    (voltage_str_infile in files_list[i][1]):
                        
            filters_39kv.append( files_list[i] )    
    
    if len(filters_39kv)==0: # se non trova misure per il filtro passa al prossimo filtro
        continue
    
    t1_39kv, t2_39kv, t3_39kv, std_39kv = plot_concentrations(filters      = filters_39kv, 
                                                              teflon_blk   = teflon_blk_39kv, 
                                                              standard     = standard_39kv, 
                                                              mylar        = mylar_39kv, 
                                                              std_conc     = std_conc_39kv, 
                                                              std_err_perc = err_perc_standard,
                                                              ax = ax, 
                                                              label = '39kv', 
                                                              min_std = min_std, 
                                                              min_sample = min_sample, 
                                                              y_units ='[$\mu$g/cm$^2$]', 
                                                              color = 'gray')
    
    t2_39kv=    t2_39kv.round({ '2-FWHM_Area_tef_blk':0,'2-FWHM_Area':0, 'conc[ug/cm^2]':2,'conc_err_log[ug/cm^2]':2,\
                               'MDL[ug/cm^2]':2, 'standard_conc':1,'2fwhm_area':0,'mylar_2fwhm_area':0})
    t2_39kv.to_csv(res_folder + filenm_prefix.replace('_STD0','') + '_'+fil+'.csv', sep=' ', columns=['Z','element', 'line', \
                                        'conc[ug/cm^2]', 'conc_err_log[ug/cm^2]', 'MDL[ug/cm^2]','2-FWHM_Area', '%Fit_Err', '2-FWHM_Area_tef_blk',\
                                        'standard_conc','2fwhm_area_std','2fwhm_area_%fit_err_std','mylar_2fwhm_area'], index=False)
    ##### store MDL values for each filter

    ## 39kV
    if mdl_frame_39kv.empty:
        mdl_frame_39kv = t2_39kv[['element','line','MDL[ug/cm^2]']]
        mdl_frame_39kv = mdl_frame_39kv.rename(columns={'MDL[ug/cm^2]':'MDL[ug/cm^2]_'+fil})
    else:   
        mdl_frame_39kv = pd.merge(mdl_frame_39kv, t2_39kv[['element','line','MDL[ug/cm^2]']], on = ['element','line'])
        mdl_frame_39kv = mdl_frame_39kv.rename(columns={'MDL[ug/cm^2]':'MDL[ug/cm^2]_'+fil})


    ############### plot PIXE ###################
    pixe_df     = pd.read_csv('./res_PIXE/' + pixe_filenm, sep =' ', na_values=pixe_nan)
    
    # leggi dataframe errori ed appendi sotto dataframe misure
    pixe_err_df = pd.read_csv('./res_PIXE/' + pixe_err_filenm, sep =' ', na_values=pixe_nan)
    pixe_err_df['Spectrum'] =  pixe_err_df['Spectrum']+'_err'
    pixe_df = pd.concat([pixe_df, pixe_err_df])
    
    # set index, trasponi ed inserisci colonne element e detector
    pixe_df.set_index('Spectrum', inplace=True)
    pixe_df = pixe_df.transpose()
    pixe_df.insert(0,'element', pixe_df.index)
    pixe_df['element'] = pixe_df['element'].apply(lambda x: x.split('_')[0])
    pixe_df.insert(1,'detector',pixe_df.index)
    pixe_df['detector'] = pixe_df['detector'].apply(lambda x: x.split('_')[1])
    pixe_df = pixe_df[['element','detector', fil.upper(), fil.upper()+'_err']]
    ###### calcola errore PIXE ########
    pixe_df.insert(len(pixe_df.columns), fil+'_err_abs', 
                   ((pixe_df[fil.upper()+'_err']*0.01*pixe_df[fil.upper()])**2 + (err_perc_standard*pixe_df[fil.upper()])**2 )**0.5 )
     
    t2_elements =  t2_39kv[t2_39kv['conc[ug/cm^2]'].notna() ]['element']
    for ele in pixe_df['element']: # remove elements measured by PIXE but not by XRF
        present=False
        for e in t2_elements:
            if e==ele:
                present = True
                
        if not present:
            pixe_df = pixe_df.drop(pixe_df[pixe_df['element']==ele].index)
            

    pixe_small = pixe_df[pixe_df['detector']=='S']
    pixe_big = pixe_df[pixe_df['detector']=='B']

    # ax.scatter(pixe_small.element, pixe_small[fil]/1000, c='#35b779' , s=30, marker='s', label='SMALL', zorder=2)
    # ax.scatter(pixe_big.element,   pixe_big[fil]/1000,   c='#440154' , s=50, marker='P', label='BIG', zorder=2)

    ax.errorbar(pixe_small.element, pixe_small[fil.upper()]/1000, yerr= pixe_small[fil+'_err_abs']/1000, capsize=7, fmt = '_', c='#35b779' , label='SMALL', zorder=2)
    ax.errorbar(pixe_big.element  , pixe_big[fil.upper()]/1000  , yerr= pixe_big[fil+'_err_abs']/1000  , capsize=7, fmt = '*', c='#440154' , label='BIG', zorder=2)

    
    ax.set_yscale('log')
    ax.legend(ncol=2)
    fig.suptitle(fil)
    fig.savefig(res_folder + filenm_prefix.replace('_STD0','') + '_conc+mdl+pixe_'+fil+'.png', format='png', dpi=300)
             
    ############### plot rapporto PIXE/XRF ###############
    
    pixe_frame = pd.concat([pixe_small, pixe_big])
    pixe_xrf_frame = pd.merge(t2_39kv, pixe_frame, on='element', how = 'outer')
    
    fig, ax = plt.subplots(1,1, figsize=(6,8))
    ax.scatter(pixe_xrf_frame.element, pixe_xrf_frame['conc[ug/cm^2]']/pixe_xrf_frame[fil.upper()] * 1000, label ='conc$_{xrf}$/conc$_{pixe}$')
    ax.grid(which='both')
    ax.set_ylim([0,11])
    ax.legend()
    # maxim = round(max(isfinite(list(pixe_xrf_frame['conc[ug/cm^2]']/pixe_xrf_frame[fil] * 1000))))+1
    ax.set_yticks(linspace(0,11,23))
    ax.set_title(fil)
    # ax.set_yscale('log')
    plt.savefig(res_folder + filenm_prefix.replace('_STD0','')+'_'+fil+'_pixe_xrf_ratio.png', format='png', dpi=300)
    plt.show()
    print('\n\n',fil, " teflon blank contribution to 2fwhm area [%]:\n",round((t2_39kv.set_index('element')['2-FWHM_Area_tef_blk']/ t2_39kv.set_index('element')['2-FWHM_Area'])*100,2))

mdl_frame_39kv.to_csv(filenm_prefix.replace('_STD0','') + '_mdl_teflon_what_39kv', sep=' ', index=False)

# t2_39kv=t2_39kv.set_index('element')
# print('\n concentrazioni:\n',t2_39kv['conc[ug/cm^2]'].astype(str)+' +- ' + t2_39kv['conc_err_log[ug/cm^2]'].round(2).astype(str))


# stampa comandi bash per rinominare spettri gupix che sono numerati 1-28
# for file in files_list:
#     print('mv', str(int(file[0].split('.')[0].split('_')[-1][-2:]))+'.jpg', file[1]+'.jpg')



















