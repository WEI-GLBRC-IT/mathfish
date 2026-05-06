import sys
import os
import re
import numpy as np

#import importlib
#import utils  # your module
#importlib.reload(utils)

from utils import *
from series import *


### Input
Nmax=50
accode=sys.argv[1]
vars=read_vars(accode)
        

try:
#if True:
    ### Preps
    Target = re.sub(' ','',vars['organism'].upper())
    Targetori = re.sub('U', 'T',Target)
    ts=[int(vars['ts1']),int(vars['ts2'])]
    L=int(vars['pl'])
    
    Nprobes=ts[1]-ts[0]-L+2
    if Nprobes>Nmax:
        message=f'ERROR: The entry triggers the evaluation of {Nprobes} probes, but the maximum number allowed is {Nmax}.'
    else:
        DG1,TS,Probe=dG1_series(vars,accode,Targetori,ts,L)
        writeAndPlot(DG1,TS,Probe,accode,1)

        # Output
        message=f'A total of {len(Probe)} probes were evaluated, and their corresponding ΔG°1 values are plotted below as a function of the 5′ end position in the target sequence.'
        
        
    write_vars([message],accode,dlm1=';') 
except:
#if False:
    import traceback

    input_file = f'userFiles/inputVars_{accode}.txt'
    log_file = f'log/dg1_inputVars_{accode}.txt'
    no_file_marker = f'log/no_file_found_{accode}'
    error_log = f'log/error_{accode}'

    if os.path.exists(input_file):
        os.system(f'cp "{input_file}" "{log_file}"')
    else:
        with open(no_file_marker, 'w') as f:
            f.write('')

    with open(error_log, 'w') as f:
        f.write(traceback.format_exc())

    os.system('python3 cleanup_log.py')
    
### Clean up
#os.system('rm *_'+accode+'.*')



