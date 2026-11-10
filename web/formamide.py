import sys
import os
import re
import numpy as np

#import importlib
#import utils  # your module
#importlib.reload(utils)

from utils import *


### Input
accode=sys.argv[1]
#accode='xxxxx'
vars=read_vars(accode)

if not len(accode)>1:
    raise Exception(f'Wrong accession code ({accode})')
        

try:
#if True:
    ### Preps
    T=float(vars['T'])
    NA=float(vars['NA'])
    Po=float(vars['Po'])
    targetmolecule=vars['targetmolecule']
    domain=vars['domain']
    bad_alignments=[]
    Targets = [re.sub(' ','',x.upper()) for x in vars['organisms'].split(';')]
    Probes = [re.sub('U','T',re.sub(' ','',x)) for x in vars['probes'].split(';')];
    # Refuse more than 5 pairs before any ClustalW/UNAFold run: the page allows 5, and the
    # plot only has 5 colours. Backstop for formamide.php's own cap (SYS-4247).
    if len(Probes) > 5 or len(Targets) > 5:
        raise ValueError(f'too many probe/target pairs ({len(Probes)}); the maximum is 5')
    stack, deltaH, deltaS = read_stacking_file('inputFiles/stacking.txt')
    
    alignments,DG1,DG2,DG3,Length=[],[],[],[],[]
    for i in range(len(Probes)):
        Target=Targets[i]
        probe=Probes[i]
    
        Targetori = re.sub('U', 'T',Target)

        # Alignment
        [targetsite, mm, pcomplement] = probealigner(Targetori, probe, accode)
        tar = Targetori[targetsite[0]-1:targetsite[1]]; 
        lines = probe_alignment_display(probe, tar, 'p', accode)
        alignments.append('<br>'.join(lines))
        
        ## dG2
        DG2.append(deltaG2(probe, T, NA, accode))        
        
        ## dG1
        DG1_this=fnearestneighbor(probe,stack,deltaH,deltaS,T,NA)
        if mm:
            DDG1 = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite, Targetori, pcomplement, accode);
            DG1_this=DG1_this+DDG1
            #print('************************************* ',i,DG1_this,DDG1,pcomplement,targetsite)
        DG1.append(DG1_this)
        
        ## dG3        
        if len(targetmolecule)==3:
            fiveprime,threeprime=newaligner(Targetori, targetsite[0], targetsite[1], targetmolecule, domain, accode);
            parse_dnd(accode,bad_alignments, 'target-'+str(i+1))
        else:
            with open('userFiles/target_'+accode+'.seq', 'w') as file:
                file.write('>target\n')
                file.write(f'{Targetori}\n')
                fiveprime,threeprime=targetsite
        DG3_this, length = deltaG3(fiveprime, threeprime, T, NA, accode);  
        DG3.append(DG3_this)
        Length.append(length)

    formamideSmM_multiple(Length, DG1, DG2, DG3, T, Po, accode)
    alignment_warning=get_alignment_warning(bad_alignments,vars['targetmolecule'],vars['domain'],'formamie%20curve%20generator')              
    write_vars(alignments+[alignment_warning],accode,dlm1=';') 
except:
#if False:
    import traceback

    input_file = f'userFiles/inputVars_{accode}.txt'
    log_file = f'log/formamide_inputVars_{accode}.txt'
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
    
### Clean up intermediate files (preserve outputVars, heplot, faplot, facurve)
os.system('rm -f userFiles/inputVars_'+accode+'.* userFiles/probe_'+accode+'.* userFiles/probe_'+accode+'-*_'+accode+'.* userFiles/target_'+accode+'.* userFiles/combined*_'+accode+'.* userFiles/targetsite_'+accode+'.* userFiles/antisense_'+accode+'.*')



