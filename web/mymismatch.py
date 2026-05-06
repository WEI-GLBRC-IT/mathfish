import sys
import os
import re

from utils import *


### Input
accode=sys.argv[1]
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

    NonTarget = re.sub(' ','',vars['mismatch'].upper())
    NonTargetori = re.sub('U', 'T',NonTarget)

    probe = re.sub('U','T',re.sub(' ','',vars['ipprobesequence']));
    
    bad_alignments=[]

    ### Alignment
    ## PM
    [targetsite_pm, mm_pm, pcomplement_pm] = probealigner(Targetori, probe, accode)
    tar_pm = Targetori[targetsite_pm[0]-1:targetsite_pm[1]]; 
    lines_pm = probe_alignment_display(probe, tar_pm, 'p', accode)
    alignment_pm='<br>'.join(lines_pm)
    #print('#############################',probe,targetsite_pm,targetsite_mm) ###

    ## MM
    [targetsite_mm, mm_mm, pcomplement_mm] = probealigner(NonTargetori, probe, accode)
    #tar_mm = Targetori[targetsite_mm[0]-1:targetsite_mm[1]]; # DOES NOT MATTER
    #lines_mm = probe_alignment_display(probe, tar_mm, 'p', accode)
    lines_mm = probe_alignment_display(probe, '_', 'p', accode)
    alignment_mm='<br>'.join(lines_mm)

    ### Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])
    Po=float(vars['Po'])
    stack, deltaH, deltaS = read_stacking_file('inputFiles/stacking.txt')

    ## dG2
    DG2=deltaG2(probe, T, NA, accode)

    ## dG1
    DG1=fnearestneighbor(probe,stack,deltaH,deltaS,T,NA)

    if mm_pm:
        DDG1_pm = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_pm, Targetori, pcomplement_pm, accode);
        DG1_pm=DG1+DDG1_pm
    else:
        DG1_pm=DG1

    if mm_mm:
        DDG1_mm = deltaG1mm(stack,deltaH,deltaS,T,NA, targetsite_mm, NonTargetori, pcomplement_mm, accode);
        DG1_mm=DG1+DDG1_mm
    else:
        DG1_mm=DG1    
        
    #print(DDG1_mm)    
    ## dG3
    targetmolecule,domain=vars['targetmolecule'],vars['domain']
    if len(targetmolecule)==3:
        fiveprime_pm,threeprime_pm=newaligner(Targetori, targetsite_pm[0], targetsite_pm[1], targetmolecule, domain, accode);
        parse_dnd(accode,bad_alignments,'target')
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{Targetori}\n')
            fiveprime_pm,threeprime_pm=targetsite_pm
    print('######################',targetsite_pm,fiveprime_pm,threeprime_pm)
    DG3_pm, length_pm = deltaG3(fiveprime_pm, threeprime_pm, T, NA, accode);

    if len(targetmolecule)==3:
        fiveprime_mm,threeprime_mm=newaligner(NonTargetori, targetsite_mm[0], targetsite_mm[1], targetmolecule, domain, accode);
        parse_dnd(accode,bad_alignments,'non-target')
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{NonTargetori}\n')
            fiveprime_mm,threeprime_mm=targetsite_mm
    print('######################',targetsite_mm,fiveprime_mm,threeprime_mm)
    DG3_mm, length_mm = deltaG3(fiveprime_mm, threeprime_mm, T, NA, accode);

    ## dGov
    DGov_pm=deltaGov(T, DG1_pm, DG2, DG3_pm, Po)
    he_pm=getHE(DGov_pm, Po, T)

    DGov_mm=deltaGov(T, DG1_mm, DG2, DG3_mm, Po)
    he_mm=getHE(DGov_mm, Po, T)

    ## Formamide denaturation
    fam_pm,fam_mm= formamideSmM_dual((length_pm,length_mm),(DG1_pm,DG1_mm), DG2, (DG3_pm,DG3_mm), T, Po, accode)

    ## Warning
    alignment_warning=get_alignment_warning(bad_alignments,vars['targetmolecule'],vars['domain'],'mismatch%20analysis',max([DG2,DG3_pm,DG3_mm]))
        
    #print(DGov,he)
    #print(alignment_warning)
    write_vars([alignment_pm,mm_pm,round(DG1_pm,1),round(DG2,1),round(DG3_pm,1),round(DGov_pm,1),round(he_pm,4),fam_pm]+\
               [alignment_mm,mm_mm,round(DG1_mm,1),round(DG2,1),round(DG3_mm,1),round(DGov_mm,1),round(he_mm,4),fam_mm]+\
               [round(DG1_mm-DG1_pm,2),round(DG3_mm-DG3_pm,2),round(DGov_mm-DGov_pm,2),round(he_mm-he_pm,4),compare_values(fam_pm,fam_mm)]+\
               [alignment_warning],accode,dlm1=';') #duplicated dG2 just for readability and convenient indexing
except:
#if False:
    import traceback

    input_file = f'userFiles/inputVars_{accode}.txt'
    log_file = f'log/mismatch_inputVars_{accode}.txt'
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



