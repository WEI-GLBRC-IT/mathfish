import sys
import os
import re

from utils import *


### Input
accode=sys.argv[1]
vars=read_vars(accode)

try:
    ### Preps
    try:
        targetsite=[int(vars['ts1']),int(vars['ts2'])]
    except:
        targetsite=None
    Target = re.sub(' ','',vars['organism'].upper())
    Targetori = re.sub('U', 'T',Target)
    if vars['ipprobesequence']=='NA':
        probe=None
    else:
        probe = re.sub('U','T',re.sub(' ','',vars['ipprobesequence']));

    if vars['iptargetsequence']=='NA':
        targetseq=None
    else:
        targetseq = re.sub('U','T',re.sub(' ','',vars['iptargetsequence']));

    #print("#############################",probe,targetsite) ###

    ### Alignment
    if vars['entry']=='p':
        [targetsite, mm, pcomplement] = probealigner(Targetori, probe, accode)
    elif vars['entry']=='t':
        probe=fprobeseq(Targetori,int(vars['ts1']),int(vars['ts2']));
        mm=0
    elif vars['entry']=='s':
        probe=fprobeseq(targetseq,1,len(targetseq));
        [targetsite, mm, pcomplement] = probealigner(Targetori, probe, accode);
    #print('#############################',probe,targetsite) ###
    tar = Targetori[targetsite[0]-1:targetsite[1]]; ###converted targetsite[0] to targetsite[0]-1 on 11/10/2024
    #print('#############################',tar) ###
    lines = probe_alignment_display(probe, tar, vars['entry'], accode)
    alignment='<br>'.join(lines)

    ### Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])
    Po=float(vars['Po'])
    stack, deltaH, deltaS = read_stacking_file('inputFiles/stacking.txt')

    ## dG2
    DG2=deltaG2(probe, T, NA, accode)

    ## dG1
    DG1=fnearestneighbor(probe,stack,deltaH,deltaS,T,NA)
    if mm:
        DDG1 = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite, Targetori, pcomplement, accode);
        DG1=DG1+DDG1
    #print(DDG1)    
    ## dG3
    targetmolecule,domain=vars['targetmolecule'],vars['domain']
    bad_alignments=[]
    if len(targetmolecule)==3:
        fiveprime,threeprime=newaligner(Targetori, targetsite[0], targetsite[1], targetmolecule, domain, accode);
        parse_dnd(accode,bad_alignments, 'target')
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{Targetori}\n')
            fiveprime,threeprime=targetsite
    print('######################',targetsite,fiveprime,threeprime)
    DG3, length = deltaG3(fiveprime, threeprime, T, NA, accode);

    ## dGov
    #print(DG1,DG2,DG3) #####
    DGov=deltaGov(T, DG1, DG2, DG3, Po)
    he=dGovHE(DGov, Po, T, accode)

    ## Formamide denaturation
    beta,fa,y, fam= formamideSmM(length,DG1, DG2, DG3, T, Po, accode)


    #print(DGov,he)
    #print('test5')
    alignment_warning=get_alignment_warning(bad_alignments,vars['targetmolecule'],vars['domain'],'General%20Analysis',max(DG2,DG3)) 
    write_vars([alignment,mm,round(DG1,1),round(DG2,1),round(DG3,1),round(DGov,1),round(he,4),fam,alignment_warning],accode,dlm1=';')
except:
#if False:
    import traceback

    input_file = f'userFiles/inputVars_{accode}.txt'
    log_file = f'log/general_inputVars_{accode}.txt'
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



