#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:19:55 2023

@author: cosimo
"""

from pandas import read_csv
from io import StringIO
from numpy import nan

def get_gupix_tables(filenm):
    # header riportati in output gupix
    header_1 = '# ID FINAL_VALUE INIT_VALUE ERROR %ERROR LAST_CHANGE %CHANGE'
    header_2 = '# Z ID DF_n/m Peak_Area 2-FWHM_Area %Fit_Err %Stat_Err +1%_Overlap LOD_Area'
    header_3 = 'Z ID # Area_Counts Yield/uC/ng/cm2 Det_Eff(-3) Filter_Trans(-5) Area_counts %Stat_Err %Fit_Err LOD_Counts X'
    try:
        infile = open(filenm, 'r')
    except:
        print('ERROR: no file found')
        raise Exception
    line = infile.readline()
    while ' #########' not in line:
        line = infile.readline()
    
    for i in range(3):
        infile.readline()
    
    table_1 = ''
    line    = ''
    while '----------------------' not in line:

        newline =  line[0:4] + line[4:12].replace(' ','')+line[12:] # remove white spaces in ID column

        try: # check if first column is float. If True, insert a nul value xxx/x in the first column
            float(newline.split()[1])
            aa = newline.split()
            aa.insert(1,'xxx/x')
            newline = ' '.join(aa) +'\n'
    
        except:
            pass
            
        table_1 = table_1 + newline
        line = infile.readline()

    table_1 = read_csv(StringIO(table_1), sep=' +', engine='python',
                          names=header_1.split())
    
    for i in range(3):
        infile.readline()
    
    table_2 = ''
    line    = ''
    while '----------------------' not in line:
        newline =  line[0:8] + line[8:15].replace(' ','') + ' ' + line[17:21].replace(' ','') + line[21:] # remove white spaces in ID column

        table_2 = table_2 + newline
        line = infile.readline()
    
    table_2 = read_csv(StringIO(table_2), sep=' +', engine='python',
                           names=header_2.split())

    while ' File:' not in line:
        line = infile.readline()
        
    for i in range(8):
        infile.readline()
    
    table_3 = ''
    line    = ''
    while ('----------------------' not in line) & ('(A "*" by the LOD value indicates' not in line):
        newline =  line[0:4] + line[4:7].replace(' ','')+line[7:] # remove white spaces in ID column
        if len(newline)>0:
            if newline[5] == ' ':
                lista=(list(newline))
                lista[5]=''
                newline = ''.join(lista)
            
        table_3 = table_3 + newline
        line = infile.readline()
     
    # print('\n\n\n',filenm,'\n\n\n')
    # print('\n\n\n',\table_3,'\n\n\n')
   
    table_3 = read_csv(StringIO(table_3), sep=' +', engine='python',
                           names=header_3.split(), index_col=False)
    # print('\n\n\n',table_3,'\n\n\n')
                            
    del table_1['#'], table_2['#']
    
    # insert column specifying the line (K, LA, MA). Insert column with element name
    
    for t in [table_1, table_2, table_3]:
        
        t.insert(1, 'line', nan)
        t.insert(1,'element', t['ID'])
        # fill the 'line' and 'element' column
        for l in ['K','LA','MA']:
            t.loc[t['ID'].apply(lambda x: x.split('/')[0].endswith(l)), 'line'] = l
            t['element'] = t['element'].apply(lambda x: x.split('/')[0].removesuffix(l))
        
    return table_1, table_2, table_3























