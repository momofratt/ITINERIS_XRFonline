#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 12:06:08 2023

@author: cosimo
"""
from read_gupix_output_2 import get_gupix_tables
import pandas as pd
from numpy import nan







def plot_concentrations(filters, teflon_blk, standard, mylar, std_conc, std_err_perc, ax, label, min_std, min_sample, y_units, color, marker='o'):
    
    std_2fwhm = get_std_2fwhm(standard, mylar, std_conc)
    T1        , T2        , T3         = get_gupix_tables('./results/'+filters[0][0])
    
    T2 = pd.merge(T2, std_2fwhm[['element','line','standard_conc','2fwhm_area', '2fwhm_area_%fit_err','mylar_2fwhm_area']], on=['element','line'])
    T2 = T2.rename(columns={'2fwhm_area':'2fwhm_area_std', '2fwhm_area_%fit_err':'2fwhm_area_%fit_err_std'})

    
    if len(teflon_blk)>0:
        T1_tef_blk, T2_tef_blk, T3_tef_blk = get_gupix_tables('./results/'+teflon_blk[0][0])
        T2 = pd.merge(T2,T2_tef_blk[['element','line','2-FWHM_Area']], on=['element','line'], suffixes=('','_tef_blk'))
        
    else: # se non trova misure teflon aggiunge colonna di zeri e stampa un warning
        T2.insert(len(T2.columns), '2-FWHM_Area_tef_blk', 0)
        print('\n ***************** WARNING: teflon blank file not found *****************\n',[filters[i][1] for i in range(len(filters))],'\n')

    T2.insert(len(T2.columns), 'conc[ug/cm^2]', nan)
    T2.insert(len(T2.columns), 'MDL[ug/cm^2]',  nan)
    T2['ID']=T2['ID'].apply(lambda x: x.replace('/0',''))
    
    # for element in std_2fwhm['element']:
        # r =   std_2fwhm[std_2fwhm['element']==element]['standard_conc'].iat[0] \
        #     /(std_2fwhm[std_2fwhm['element']==element]['2fwhm_area'   ].iat[0] - std_2fwhm[std_2fwhm['element']==element]['mylar_2fwhm_area'].iat[0] )\
        #         * min_std / min_sample
        # T2.loc[T2['element']==element, 
        #             'conc[ug/cm^2]'] =  (T2['2-FWHM_Area']-T2['2-FWHM_Area_tef_blk'])*r
        # T2.loc[T2['element']==element, 
        #             'MDL[ug/cm^2]']  =  T2['LOD_Area']*r
        
    ################### Calcolo errore con derivate parziali (da rivedere!!!!) ################################
    # D  = T2['2fwhm_area_std'] 
    # M =  min_std / min_sample
    # A1 =  T2['%Fit_Err']    *0.01           *T2['2-FWHM_Area'] * T2['standard_conc'] * M / D
    # A2 =  T2['2fwhm_area_%fit_err_std']*0.01*T2['2fwhm_area_std']*T2['2-FWHM_Area'] * T2['standard_conc'] * M  / (D**2)
    # A3 =  T2['2-FWHM_Area'] * std_err_perc*T2['conc[ug/cm^2]'] * T2['standard_conc'] * M / D
    # T2.insert(len(T2.columns), 'A1', A1)
    # T2.insert(len(T2.columns), 'A2', A2)
    # T2.insert(len(T2.columns), 'A3', A3)
    # conc_error = A1 + A2 + A3
    # T2.insert(T2.columns.get_loc('conc[ug/cm^2]')+1, 'conc_err[ug/cm^2]' , conc_error)
    ###########################################################################################################
    
    r = T2['standard_conc'] / (T2['2fwhm_area_std'] - T2['mylar_2fwhm_area'] ) * min_std / min_sample
    T2['conc[ug/cm^2]'] =  (T2['2-FWHM_Area']-T2['2-FWHM_Area_tef_blk'])*r
    T2['MDL[ug/cm^2]']  =   T2['LOD_Area']*r
        
    T2 = T2[T2['conc[ug/cm^2]'].notna()]
    T2 = T2[T2['conc[ug/cm^2]']>=0]    
    
    ################## Calcolo errore con derivate log ##############################
    err_log = (T2['%Fit_Err']                * 0.01 + \
               T2['2fwhm_area_%fit_err_std'] * 0.01 + \
               std_err_perc) * T2['conc[ug/cm^2]'] 
    
    T2.insert(T2.columns.get_loc('conc[ug/cm^2]')+1, 'conc_err_log[ug/cm^2]' , err_log)
    
    # ax.errorbar(T2['element'],
    #             T2['conc[ug/cm^2]'], 
    #             yerr=( (T2['conc[ug/cm^2]']*T2['%Fit_Err']*0.01)**2 + (T2['conc[ug/cm^2]']*std_err_perc)**2)**0.5, 
    #             capsize=5, fmt = marker, label=label, color=color, zorder=3)
   
    ax.errorbar(T2['element'],
                T2['conc[ug/cm^2]'], 
                yerr = T2['conc_err_log[ug/cm^2]'], 
                capsize=5, fmt = marker, label=label, color=color, zorder=3)
    
    ax.scatter(T2['element'],
                T2['MDL[ug/cm^2]'], 
                marker = '_', color=color, zorder=3)
        
    ax.grid(True ,which='both', zorder=0)
    ax.legend(ncol=3)
    ax.set_ylabel('Conc '+y_units)
    for label in ax.get_xticklabels(): # format last axis xticks
        label.set_ha("center")
        label.set_rotation(45)
    
    return T1,T2,T3, std_2fwhm

def get_std_2fwhm(std_list, mylar, std_conc, area = '2-FWHM_Area'):
    """ ritorna lista con le misure di area degli standard data una lista di files. 
    Lo standard è la 2-FWHM_Area ma può essere richiesta anche la peak area """
    
    area_str = area.replace('-','').lower()
    
    std_2fwhm = pd.DataFrame(data={'element':[], 
                                   'line':[], 
                                   area_str:[],
                                   area_str+'_%fit_err':[],
                                   'LOD':[],
                                   'mylar_'+area_str:[], 
                                   'mylar_lod':[], 
                                   'standard_conc':[], 
                                   'Z':[]})
    
    _, t2_myl, _ = get_gupix_tables('./results/'+mylar[0][0])
   
    for s_file in std_list:
        
        element = s_file[1].split('_')[1]
        if element.startswith('K') | (element=='P') | element.startswith('V'):
            element=element[0]
        elif element.startswith('Mylar'):
            element=element
        elif element == 'CuSx':
            element='S'
        else:
            element=element[0:2]

        line = std_conc[std_conc['ele' ]==element]['line'].iat[0]

        _, t2, _ = get_gupix_tables('./results/'+s_file[0])
        # appendi riga al dataframe per ogni elemento

        condition = (t2['line']==line) & (t2['element'] == element) 
        condition_myl = (t2_myl['line']==line) & (t2_myl['element'] == element) 

        # inserisce riga nel dataframe std_2fwhm
        std_2fwhm.loc[len(std_2fwhm.index)] = [element, line,
                            t2[                      condition][area].iat[0],
                            t2[                      condition]['%Fit_Err'].iat[0],
                            t2[                      condition]['LOD_Area'].iat[0],
                            t2_myl[              condition_myl][area].iat[0],
                            t2_myl[              condition_myl]['LOD_Area'].iat[0],
                            std_conc[std_conc['ele' ]==element]['conc'].iat[0],
                            t2[                      condition]['Z'].iat[0]]
            
    return std_2fwhm

def merge_and_plot(t2, t3, ax, filenm, legend, color, marker, ylab_str='', 
                   text=False, LOD=True, LOD_color=False, custom_lab=''):
    """
    Perform merge with t3 and plot on a given axis

    Parameters
    ----------
    t2 : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.
    filenm : TYPE
        DESCRIPTION.
    legend : TYPE
        DESCRIPTION.
    color : TYPE
        DESCRIPTION.
    marker : TYPE
        DESCRIPTION.
    LOD: bool
        whether or not to plot LOD
    LOD_color:
        whether or not to plot coloured or black LOD

    Returns
    -------
    ax : TYPE
        DESCRIPTION.
    t2 : TYPE
        DESCRIPTION.

    """
    
    t2 = t2.merge(t3[['element','line','X']], on=['element','line'], how='left')
    t2.sort_values('Z')
    t2.insert(len(t2.columns),'?', t2['2-FWHM_Area'])
    t2.loc[ t2['X'] != '?', '?' ] = nan
    t2.loc[ t2['X'] != 'Y', '2-FWHM_Area' ] = nan
    if custom_lab == '':
        label = ' '.join(filenm.split('_')[0:2])
    else:
        label = custom_lab
        
    ax.xaxis.set_tick_params(labelbottom=True)
    
    ax.errorbar(t2['element'] + ' '+ t2['line'],
                t2['2-FWHM_Area'], 
                yerr=t2['2-FWHM_Area']*t2['%Fit_Err']*0.01, 
                label=label,
                capsize=5, fmt = 'o',
                c=color)
    
    ax.errorbar(t2['element'] + ' '+ t2['line'],
                t2['?'], 
                yerr=t2['?']*t2['%Fit_Err']*0.01, 
                label= label + ' ?',
                capsize=5, fmt = 'o', c='dark'+color)
    if LOD:
        if LOD_color:
            LOD_col = color
        else:
            LOD_col = 'black'
            
        ax.scatter(t2['element'] + ' '+ t2['line'],
                    t2['LOD_Area'], 
                    label=label + ' LOD',
                    marker = marker, c=LOD_col)
        
    # ax.set_xticks(t2['element'] + ' ' + t2['line'])
    ax.grid(visible=True,which='both')
    if legend:
        ax.legend(ncol=3, bbox_to_anchor=(0.8, 1.2), fancybox=True,shadow=True)
    ax.set_yscale('log')

    ax.set_ylabel('2.FWHM Area'+ylab_str)
    if text:
        t=ax.text(0.85,  0.05,  ' '.join(filenm.split('_')[2:]) , horizontalalignment='center', size='large', color='black',transform=ax.transAxes)                                              
        t.set_bbox(dict(facecolor='grey', alpha=0.3))
 
    for label in ax.get_xticklabels(): # format last axis xticks
        label.set_ha("center")
        label.set_rotation(45)

    return ax, t2
