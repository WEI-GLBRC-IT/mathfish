import sys
import os
import re
import numpy as np

import importlib
import utils  # your module
importlib.reload(utils)

from utils import *


### Input
accode=sys.argv[1]
#accode='aziz'
vars=read_vars(accode)

### Functions
def compare_values(var1, var2):
    def parse(val):
        if isinstance(val, (int, float)):
            return val
        elif isinstance(val, str):
            if val.startswith("<"):
                return ("<", float(val[2:]))
            elif val.startswith(">"):
                return (">", float(val[2:]))
        return None  # unexpected format

    # Case 1: Both are numbers
    if isinstance(var1, (int, float)) and isinstance(var2, (int, float)):
        return round(var2 - var1, 2)

    # Case 2: Both are the same string
    if isinstance(var1, str) and isinstance(var2, str) and var1 == var2:
        return 'NaN'

    # Parse values
    p1 = parse(var1)
    p2 = parse(var2)

    # Case 3: One or both are strings
    try:
        if isinstance(var1, (int, float)):
            if p2[0] == "<":
                return f"< {round(p2[1] - var1, 2)}"
            elif p2[0] == ">":
                return f"> {round(p2[1] - var1, 2)}"
        elif isinstance(var2, (int, float)):
            if p1[0] == "<":
                return f"< {round(var2 - p1[1], 2)}"
            elif p1[0] == ">":
                return f"< {round(var2 - p1[1], 2)}"
        else:
            # both are strings but not equal
            if p1[0] == ">" and p2[0] == "<":
                return f"< {round(p2[1] - p1[1], 2)}"
            elif p1[0] == "<" and p2[0] == ">":
                return f"> {round(p2[1] - p1[1], 2)}"
            else:
                return 'NaN'
    except Exception:
        return 'NaN'
        


try:
#if True:
    ### Preps
    Target = re.sub(' ','',vars['organism'].upper())
    Targetori = re.sub('U', 'T',Target)

    NonTarget = re.sub(' ','',vars['nontarget'].upper())
    NonTargetori = re.sub('U', 'T',NonTarget)

    probe = re.sub('U','T',re.sub(' ','',vars['probe']));
    competitor = re.sub('U','T',re.sub(' ','',vars['compe']));

    ## Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])
    Po=float(vars['Po'])
    Co=float(vars['Co'])
    stack, deltaH, deltaS = read_stacking_file('inputFiles/stacking.txt')    
    
    bad_alignments=[]

    #### PROBE (alignment,DG1,DG2)
    ### Alignment
    ## Probe-target
    [targetsite_pt, mm_pt, pcomplement_pt] = probealigner(Targetori, probe, accode)
    #tar_pt = Targetori[targetsite_pt[0]-1:targetsite_pt[1]]; 
    lines_pt = probe_alignment_display(probe, '_', 'p', accode)
    alignment_pt='<br>'.join(lines_pt)
    #print('#####PT\n'+'\n'.join(lines_pt)+'/n#####') ###

    ## Probe-nontarget
    [targetsite_pn, mm_pn, pcomplement_pn] = probealigner(NonTargetori, probe, accode)
    #tar_pn = Targetori[targetsite_pn[0]-1:targetsite_pn[1]]; 
    lines_pn = probe_alignment_display(probe, '_', 'p', accode)
    alignment_pn='<br>'.join(lines_pn)
    #print('#####PN\n'+'\n'.join(lines_pn)+'/n#####') ###

    ## dG2
    DG2_p=deltaG2(probe, T, NA, accode)

    ## dG1
    DG1_p=fnearestneighbor(probe,stack,deltaH,deltaS,T,NA)
    
    if mm_pt:
        DDG1_pt = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_pt, Targetori, pcomplement_pt, accode);
        DG1_pt=DG1_p+DDG1_pt
    else:
        DG1_pt=DG1_p

    if mm_pn:
        DDG1_pn = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_pn, NonTargetori, pcomplement_pn, accode);
        DG1_pn=DG1_p+DDG1_pn
        print('here',DDG1_pn)
    else:
        DG1_pn=DG1_p  


    #### COMPETITOR (alignment,DG1,DG2)
    ### Alignment
    ## Competitor-target
    [targetsite_ct, mm_ct, ccomplement_ct] = probealigner(Targetori, competitor, accode)
    #tar_ct = Targetori[targetsite_pt[0]-1:targetsite_pt[1]]; 
    lines_ct = probe_alignment_display(competitor, '_', 'p', accode)
    alignment_ct='<br>'.join(lines_ct)
    #print('#####CT\n'+'\n'.join(lines_ct)+'/n#####') ###

    ## Competitor-nontarget
    [targetsite_cn, mm_cn, ccomplement_cn] = probealigner(NonTargetori, competitor, accode)
    #tar_pn = Targetori[targetsite_pn[0]-1:targetsite_pn[1]]; 
    lines_cn = probe_alignment_display(competitor, '_', 'p', accode)
    alignment_cn='<br>'.join(lines_cn)
    #print('#####CN\n'+'\n'.join(lines_cn)+'/n#####') ###

    ## dG2
    DG2_c=deltaG2(competitor, T, NA, accode)
    
    DG2=[DG2_p,DG2_p,DG2_c,DG2_c]
    #print(DG2)    
    
    ## dG1
    DG1_c=fnearestneighbor(competitor,stack,deltaH,deltaS,T,NA)
    
    if mm_ct:
        DDG1_ct = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_ct, Targetori, ccomplement_ct, accode);
        DG1_ct=DG1_c+DDG1_ct
    else:
        DG1_ct=DG1_c

    if mm_cn:
        DDG1_cn = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_cn, NonTargetori, ccomplement_cn, accode);
        DG1_cn=DG1_c+DDG1_cn
    else:
        DG1_cn=DG1_c 
    
    DG1=[DG1_pt,DG1_pn,DG1_ct,DG1_cn]
    #print(DG1)
        
    #print(DDG1_mm)    
    
    #### dG3
    ### Target    
    targetmolecule,domain=vars['targetmolecule'],vars['domain']
    if len(targetmolecule)==3:
        fiveprime_pt,threeprime_pt=newaligner(Targetori, targetsite_pt[0], targetsite_pt[1], targetmolecule, domain, accode);
        fiveprime_ct,threeprime_ct=newaligner(Targetori, targetsite_ct[0], targetsite_ct[1], targetmolecule, domain, accode,skip_alignment=True);
        parse_dnd(accode,bad_alignments,'target')
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{Targetori}\n')
            fiveprime_pt,threeprime_pt=targetsite_pt
            fiveprime_ct,threeprime_ct=targetsite_ct
    print('######################',targetsite_pt,fiveprime_pt,threeprime_pt)
    print('######################',targetsite_ct,fiveprime_ct,threeprime_ct)
    DG3_pt, length_pt = deltaG3(fiveprime_pt, threeprime_pt, T, NA, accode);
    DG3_ct, length_ct = deltaG3(fiveprime_ct, threeprime_ct, T, NA, accode);

    ### NonTarget
    if len(targetmolecule)==3:
        fiveprime_pn,threeprime_pn=newaligner(NonTargetori, targetsite_pn[0], targetsite_pn[1], targetmolecule, domain, accode);
        fiveprime_cn,threeprime_cn=newaligner(NonTargetori, targetsite_cn[0], targetsite_cn[1], targetmolecule, domain, accode,skip_alignment=True);
        parse_dnd(accode,bad_alignments,'non-target')
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{NonTargetori}\n')
            fiveprime_pn,threeprime_pn=targetsite_pn
            fiveprime_cn,threeprime_cn=targetsite_cn
    print('######################',targetsite_pn,fiveprime_pn,threeprime_pn)
    print('######################',targetsite_cn,fiveprime_cn,threeprime_cn)
    DG3_pn, length_pn = deltaG3(fiveprime_pn, threeprime_pn, T, NA, accode);
    DG3_cn, length_cn = deltaG3(fiveprime_cn, threeprime_cn, T, NA, accode);

    DG3=[DG3_pt,DG3_pn,DG3_ct,DG3_cn]
    length=[length_pt,length_pn,length_ct,length_cn]
    #print(DG3)

    #### dGov
    DGov=4*[np.nan]
    Oo=[Po,Po,Co,Co]
    for i in (0,1,2,3):
        DGov[i]=deltaGov(T, DG1[i], DG2[i], DG3[i], Oo[i])
        
    #### Hybridization Efficiency
    he=4*[np.nan]
    for i in (0,1,2,3):
        he[i]=getHE(DGov[i], Oo[i], T)


    #print(DGov)
    #print(he)

    ## Formamide denaturation
    formamideSmM_competitor(length,DG1, DG2, DG3, T, Oo, accode)

    ## Warning
    alignment_warning=get_alignment_warning(bad_alignments,vars['targetmolecule'],vars['domain'],'competitor%20analysis',max(DG2+DG3))
        
    #print(DGov,he)
    #print(alignment_warning)
    write_vars([alignment_pt,alignment_pn,alignment_ct,alignment_cn]+[round(x,1) for x in DG1]+[round(x,1) for x in DG2]+\
               [round(x,1) for x in DG3]+[round(x,1) for x in DGov]+[round(x,2) for x in he]+[alignment_warning],accode,dlm1=';') 
except:
#if False:
    import traceback

    input_file = f'userFiles/inputVars_{accode}.txt'
    log_file = f'log/compe_inputVars_{accode}.txt'
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



