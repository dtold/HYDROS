#!/usr/bin/env python
from __future__ import print_function, division
import sys
import os
import numpy as np
np.set_printoptions(threshold='nan')

#add source code directory to path
startdir=os.path.join(os.path.split(sys.argv[0])[0])
sys.path+=[os.path.join(startdir,'src')]

from scan_tools import HYDROS

#read input parameters from file name given on command line
try: 
    inputfile=sys.argv[1]
except:
    inputfile='parameters'

par=open(inputfile,'r')
p={}
for line in par.readlines():
    line=line.strip()
    if len(line)>0:
        if line[0]!='#':
            if '#' in line: line=line[:line.index('#')]
            var,val=line.split('=')
            p[var]=val.strip()

#from M. Oberparleitner's additions to ParIO.py in GENE:
for var in p:
    try:   
        p[var] = int(p[var])
    except ValueError:
        try:  
            p[var] = float(p[var])
        except ValueError:
            pass
for var in p:
    try:
        p[var]=eval(p[var])
    except:
        pass

#set parameters as used in the code
params=[p['wavevector_mode'],p['kperp'],p['kpar'],p['beta'],p['tau'],p['Tpar_Tperp'],p['gam'],p['eta'],p['nb'],p['theta'],p['k']]

#here we prepare and run the actual scan
scan_options= [p['radius'],p['dev_limit'],p['cp_limit']]
scan1=        [p['scanvar'],p['nval'],p['scan_range'],p['log'],p['scan_list1']]
scan2=        [p['scanvar2'],p['nval2'],p['scan_range2'],p['log2'],p['scan_list2']]

#generate output file paths
resultfilepath=os.path.join('output/',p['resultfile'])
logfilepath=os.path.join('output/',p['logfile'])

#run the code
HYDROS(params,p['start'],scan1,scan2,p['log'],p['method'],scan_options,
       scan_type=p['scan_type'],outfilename=logfilepath,scanfilename=resultfilepath,
       normalization=p['normalization'])


print('_'*80+'\n')
print('Output files: %s, %s\n\n' %(resultfilepath,logfilepath))
