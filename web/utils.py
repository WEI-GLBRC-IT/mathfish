import sys
import os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import math

R=1.99E-3 #Gas constant
maxDist=0.25 #Maximum tolerable distance between aligned sequences (for dG3 calculations) before a warning is issued.
ref_organisms={'a':('Archaea','Methanosarcina barkeri'),'b':('Bacteria','Escherichia coli'),'e':('Eukarya','Saccharomyces cerevisiae')}
max_folding_dG=5. #Maximum allowed folding free energy (deltaG2 or deltaG3)

def read_vars(access_code):
    '''Reads input variables provided in a file saved by PHP'''
    filepath = 'userFiles/inputVars_' + access_code + '.txt'
    try:
        with open(filepath, 'r') as file:
            variables = {}
            current_var_name = None
            content_lines = []
            for line in file:
                if 'VAR_' in line:
                    if current_var_name:
                        variables[current_var_name] = ''.join(content_lines).strip()
                        content_lines = []
                    current_var_name = line.strip().split('_')[1]
                else:
                    content_lines.append(line.strip())
            if current_var_name:
                variables[current_var_name] = ''.join(content_lines).strip()
    except FileNotFoundError:
        print(f"Error: The file {filepath} does not exist.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    return variables

def write_vars(vars,access_code,dlm1=';',dlm2=','):
    '''Writes variables in list <vars> to file outputVars_access_code as separated by <dlm1>. 
    If a variable is a list, then the elements will further be separated by <dlm2>.
    '''
    with open('userFiles/outputVars_'+access_code+'.txt', 'w') as file:
        for var in vars:
            if isinstance(var,list):
                var=dlm2.join([str(x) for x in var])
            file.write(str(var))
            file.write(dlm1)

def read_stacking_file(filename):
    # Initialize lists to store the data
    stack = []
    deltaH = []
    deltaS = []

    # Open the file and read the contents
    with open(filename, 'r') as file:
        for line in file:
            # Split the line into components based on tabs
            parts = line.strip().split('\t')
            # Append data to respective lists
            stack.append(parts[0])
            deltaH.append(float(parts[1]))
            deltaS.append(float(parts[2]))

    return stack, deltaH, deltaS

def get_DG_from_file(fname):
    # Get the parameters from UNAFold output
    with open(fname, 'r') as file:
        lines = file.readlines()
        # Skip the header line and extract the second line for delta G value
        DG = float(lines[1].split()[1])  # Assuming the delta G value is the second item on the second line
    return DG
    
def get_lines_from_file(filename):
    # Read lines from a file without stripping whitespace
    with open(filename, 'r') as file:
        lines = file.readlines()  # This reads all lines with trailing newlines
    return lines
    
    
def fprobeseq(target, fiveprime, threeprime):
    # Convert target to uppercase
    target = target.upper()
    
    # Adjust the indices for Python's 0-based indexing and handle reverse slicing
    if fiveprime == 1:
        template = target[threeprime-1::-1]  # Go from threeprime to the start
    else:
        template = target[threeprime-1:fiveprime-2:-1]  # Normal reverse slice ###converted fiveprime-1 to fiveprime-2 on 11/10/2024
    
    # Define the DNA matching pairs
    dna = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Create the probe by replacing each nucleotide with its complement
    probe = ''.join(dna[nuc] for nuc in template)

    return probe
    
def cumsum_non_gaps(aligned_seq):
    '''Subroutine of probealigner'''
    cumsum = []
    count = 0
    for char in aligned_seq:
        if char != '-':
            count += 1
        cumsum.append(count)
    return cumsum

def probealigner(target_ori, probe, access_code):
    # Convert probe to its complement
    probe = probe.upper()
    p_complement = fprobeseq(probe, 1, len(probe))

    # Combine the sequences in one file for alignment
    with open('userFiles/combinedprobetarget_'+access_code+'.txt', 'w') as f:
        f.write('>T\n')
        f.write(f'{target_ori}\n')
        f.write('>P\n')
        f.write(f'{p_complement}\n')

    # Call the alignment software
    os.system('clustalw2 userFiles/combinedprobetarget_'+access_code+'.txt /output=GDE /case=upper')

    # Read from the GDE formatted output file
    with open('userFiles/combinedprobetarget_'+access_code+'.gde', 'r') as f:
        lines = f.readlines()

    # Parse the alignment results
    target_aligned = ''
    probe_aligned = ''
    in_target = False
    in_probe = False
    for line in lines:
        if line.startswith('#'):
            in_target = 'T' in line
            in_probe = 'P' in line
            continue
        if in_target:
            target_aligned += line.strip()
        elif in_probe:
            probe_aligned += line.strip()
    
    # Determining target site and mismatches
    target_net = cumsum_non_gaps(target_aligned)
    probe_net = cumsum_non_gaps(probe_aligned)

    # Finding the starting and ending match positions
    if probe_net:
        start_index = next((i for i, x in enumerate(probe_net) if x == 1), None)
        end_index = next((i for i, x in enumerate(probe_net) if x == len(probe)), None)
        target_site = (target_net[start_index], target_net[end_index])

    mismatches, bulges = 0, 0
    for t, p in zip(target_aligned[start_index:end_index+1], probe_aligned[start_index:end_index+1]):
        if t != p:
            mismatches += 1
            if '-' in (t, p):
                bulges += 1

    print(f"Target site position (5' and 3' terminals): {target_site}")
    print(f"The probe has {mismatches} mismatches to the target.")
    if bulges > 0:
        bulge_text = 'are bulged mismatches' if bulges > 1 else 'is a bulged mismatch'
        print(f"Of which {bulges} {bulge_text}.")

    return target_site, mismatches, p_complement

def probe_alignment_display(probe, tar, entry, access_code):
    L = len(probe)
    proberev = probe[::-1]
    line1 = ['.' for _ in range(9)]
    line2 = list('PROBE--3\'')
    line3 = ['.' for _ in range(9)]
    line4 = list('TARGET-5\'')
    line5 = ['.' for _ in range(9)]
  
    if entry == 't':
        for i in range(L):
            line2.append(proberev[i])
            line4.append(tar[i])
            line3.append('|')
            line1.append('.')
            line5.append('.')
    else:
        with open('userFiles/combinedprobetarget_'+access_code+'.gde', 'r') as file:
            lines = file.readlines()

        # Parse the alignment results
        target_align = ''
        probe_align = ''
        in_target = False
        in_probe = False
        for line in lines:
            if line.startswith('#'):
                in_target = 'T' in line
                in_probe = 'P' in line
                continue
            if in_target:
                target_align += line.strip()
            elif in_probe:
                probe_align += line.strip()

        # Using cumulative sum to find start and end index
        probe_net = cumsum_non_gaps(probe_align)
        start_index = probe_net.index(1)  # first non-gap position
        end_index = probe_net.index(len(probe)) # last non-gap position


        target = target_align[start_index:end_index+1]
        probe = probe_align[start_index:end_index+1]
        

        # Complement of reversed probe segment
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-'}
        pcomprev = [complement.get(x, '.') for x in probe]

        for i in range(end_index - start_index + 1):
            if probe[i] == target[i]:
                line1.append('.')
                line2.append(pcomprev[i])
                line4.append(target[i])
                line3.append('|')
                line5.append('.')
            else:
                line2.append('.')
                line4.append('.')
                line3.append('.')
                line1.append(pcomprev[i])
                line5.append(target[i])

    # Join lists to strings for final display
    return ''.join(line1), ''.join(line2)+"5'", ''.join(line3), ''.join(line4)+"3'", ''.join(line5)
    
def deltaG2(probe, T, NA, access_code):
    # Create sequence file
    with open('userFiles/probe_'+access_code+'.seq', 'w') as file:
        file.write('>probe\n')
        file.write(f'{probe}\n')

    # Build and implement the UNAFold command
    com1 = f"hybrid-ss-min userFiles/probe_{access_code}.seq -n DNA -t {T} -T {T} -N {NA} -o userFiles/probe_{access_code}"
    os.system(com1)

    # Get the parameters from UNAFold output
    DG2=get_DG_from_file('userFiles/probe_'+access_code+'.dG')

    return min(DG2,max_folding_dG)
    
def fnearestneighbor(mprobe, stack, deltaH, deltaS, T, NA):
    from math import log

    # p is the number of parameters
    p = len(stack)

    # Coefficients, a vector equivalent in length to stack, indicates the multipliers
    # necessary for each stacking
    coef = np.zeros(p)
    coef[-1]=1
    
    # Calculate the coefficients for each stack entry except the last one
    index_map = {nn: idx for idx, nn in enumerate(stack)}
    for i in range(len(mprobe)-1):
        nn = mprobe[i:i+2]
        if nn in index_map:
            coef[index_map[nn]] += 1
    
    # Dot products. DS additionally has salt correction.
    DH = np.dot(coef, deltaH)
    DS = np.dot(coef, deltaS) + 0.368 * (sum(coef) - 1) * log(NA)
    
    # Final calculation of ΔG
    DG = DH - ((T + 273.15) / 1000) * DS
    
    return DG

def deltaG1mm(stack, deltaH, deltaS, T, NA, targetsite, Targetori, pcomplement,access_code):
    # Read stacking parameters for DNA and RNA from files
    stack_dna, deltaH_dna, deltaS_dna = read_stacking_file('inputFiles/stackingDNA.txt')
    stack_rna, deltaH_rna, deltaS_rna = read_stacking_file('inputFiles/stackingRNA.txt')
    print(targetsite)
    print(Targetori[targetsite[0]-1:targetsite[1]])
    
    # Generate file names - full paths for file I/O, basenames for hybrid-2s.pl
    targetsite_fname,antisense_fname,probe_fname=['userFiles/'+x+'_'+access_code for x in ('targetsite','antisense','probe')]
    probe_antisense_fname,probe_targetsite_fname=['userFiles/probe_'+access_code+'-'+x+'_'+access_code for x in ('antisense','targetsite')]

    # Basenames for hybrid-2s.pl (runs from userFiles/ directory)
    targetsite_base,antisense_base,probe_base=[x+'_'+access_code for x in ('targetsite','antisense','probe')]

    # Writing target site to a file
    with open(targetsite_fname+'.seq', 'w') as file:
        file.write('>targetsite\n')
        file.write(f"{Targetori[targetsite[0]-1:targetsite[1]]}\n")

    # Writing perfect complement to a file
    with open(antisense_fname+'.seq', 'w') as file:
        file.write('>antisense\n')
        file.write(f"{pcomplement}\n")

    ## DNA Analysis
    # Execute hybrid-2s for DNA/DNA duplex (run from userFiles/ for correct output location)
    subprocess.run(f"hybrid-2s.pl -n DNA -t {T} -T {T} {probe_base}.seq {antisense_base}.seq", shell=True, cwd="userFiles")
    DG1_DNA_DNA = get_DG_from_file(probe_antisense_fname+'.dG') #B

    # Nearest neighbor calculations for DNA
    lines = get_lines_from_file(probe_antisense_fname+'.asc')
    mprobe = lines[1]
    nn_DNA_RNA_perfect = fnearestneighbor(mprobe, stack, deltaH, deltaS, T, NA) #A
    nn_DNA_DNA_perfect = fnearestneighbor(mprobe, stack_dna, deltaH_dna, deltaS_dna, T, NA) #B'

    # Mismatch calculation for DNA
    subprocess.run(f"hybrid-2s.pl -n DNA -t {T} -T {T} {probe_base}.seq {targetsite_base}.seq", shell=True, cwd="userFiles")
    DG1_DNA_DNA_mm = get_DG_from_file(probe_targetsite_fname+'.dG') #C

    # Nearest neighbor calculations excluding mismatched positions for DNA
    lines = get_lines_from_file(probe_targetsite_fname+'.asc')
    mprobe = lines[1]
    nn_DNA_DNA_mm = fnearestneighbor(mprobe, stack_dna, deltaH_dna, deltaS_dna, T, NA) #D
    nn_DNA_RNA_mm_wrt_dna = fnearestneighbor(mprobe, stack, deltaH, deltaS, T, NA) #E

    # Calculate delta G components for DNA
    DDG1_broken_wrt_DNA = nn_DNA_RNA_perfect - nn_DNA_RNA_mm_wrt_dna #F = A-E
    DDG1_DNA_broken = nn_DNA_DNA_perfect - nn_DNA_DNA_mm #G = B'- D
    DDG1_DNA_loop = DG1_DNA_DNA_mm - DG1_DNA_DNA + DDG1_DNA_broken # H = C-B+G

    ## RNA Analysis
    # Execute hybrid-2s for RNA/RNA duplex (run from userFiles/ for correct output location)
    subprocess.run(f"hybrid-2s.pl -n RNA -t {T} -T {T} {probe_base}.seq {antisense_base}.seq", shell=True, cwd="userFiles")
    DG1_RNA_RNA = get_DG_from_file(probe_antisense_fname+'.dG') #rB

    # Nearest neighbor calculations for RNA
    lines = get_lines_from_file(probe_antisense_fname+'.asc')
    mprobe = lines[1]
    nn_DNA_RNA_perfect_wrt_rna = fnearestneighbor(mprobe, stack, deltaH, deltaS, T, NA) #rA
    nn_RNA_RNA_perfect = fnearestneighbor(mprobe, stack_rna, deltaH_rna, deltaS_rna, T, NA) #rB'

    # Mismatch calculation for RNA
    subprocess.run(f"hybrid-2s.pl -n RNA -t {T} -T {T} {probe_base}.seq {targetsite_base}.seq", shell=True, cwd="userFiles")
    DG1_RNA_RNA_mm = get_DG_from_file(probe_targetsite_fname+'.dG') #rC
    
    # Nearest neighbor calculations excluding mismatched positions for RNA
    lines = get_lines_from_file(probe_targetsite_fname+'.asc')
    mprobe = lines[1]
    nn_RNA_RNA_mm = fnearestneighbor(mprobe, stack_rna, deltaH_rna, deltaS_rna, T, NA) #rD
    nn_DNA_RNA_mm_wrt_rna = fnearestneighbor(mprobe, stack, deltaH, deltaS, T, NA) #rE
    
    # Calculate delta G components for RNA
    DDG1_broken_wrt_RNA = nn_DNA_RNA_perfect_wrt_rna - nn_DNA_RNA_mm_wrt_rna #rF=rA-rE
    DDG1_RNA_broken = nn_RNA_RNA_perfect - nn_RNA_RNA_mm #rG=rB'-rD
    DDG1_RNA_loop = DG1_RNA_RNA_mm - DG1_RNA_RNA + DDG1_RNA_broken #rH=rC-rB+rG

    # Calculate final delta G values
    DDG1_loop = np.mean([DDG1_RNA_loop, DDG1_DNA_loop])
    DDG1_broken = np.mean([DDG1_broken_wrt_RNA, DDG1_broken_wrt_DNA])
    DDG1 = -DDG1_broken + DDG1_loop
    
    return DDG1

def newaligner(Targetori, targetsite1, targetsite2, targetmolecule, domain,access_code,skip_alignment=False):
    # Setting up file names based on domain and target molecule
    ref_filename = f'inputFiles/ref-{domain}{targetmolecule}.txt'
    domain_filename = f'inputFiles/domain-{domain}{targetmolecule}.txt'
    
    # Read reference original sequences
    with open(ref_filename, 'r') as file:
        Refori = file.read().replace('\n', '')    

    # Read domain borders
    with open(domain_filename, 'r') as file:
        Refdomain = []
        for line in file:
            parts = line.strip().split()
            Refdomain.append([int(parts[0]), int(parts[1])])    
    
    if not skip_alignment:
        # Combining the sequences in one file
        with open('userFiles/combined_'+access_code+'.txt', 'w') as file:
            file.write('>R\n')
            file.write(f'{Refori}\n')
            file.write('>O\n')
            file.write(f'{Targetori}\n')

        # Execute ClustalW alignment
        os.system('clustalw2 userFiles/combined_'+access_code+'.txt /output=GDE /case=upper')
    
    # Read alignment output
    alignment=''.join([x.strip() for x in get_lines_from_file('userFiles/combined_'+access_code+'.gde')])
    indexes = [i for i, x in enumerate(alignment) if x == '#']
    if len(indexes) != 2:
        raise ValueError('Alignment file is corrupted')

    Refaln = alignment[2:indexes[1]]
    Targetaln = alignment[indexes[1]+2:]    
    
    # Determine domain borders on the target
    Refnet = np.array(cumsum_non_gaps(Refaln))
    Targetnet = cumsum_non_gaps(Targetaln)
    n = len(Refdomain)
    Targetdomain = np.zeros((n, 2), dtype=int)
    #print(Refnet)
    #print(Refdomain)
    for i in range(n):
        #print(i)
        start_index = np.where(Refnet == Refdomain[i][0])[0][0]
        end_index = np.where(Refnet == Refdomain[i][1])[0][0]
        Targetdomain[i, :] = [Targetnet[start_index], Targetnet[end_index]]

    # Correction for end points
    Targetdomain[0, 0] = 1
    Targetdomain[-1, 1] = len(Targetori)
    
   # Determining the targeted domain
    indexa = max([i for i, x in enumerate(Targetdomain[:, 0]) if x <= targetsite1])
    indexb = min([i for i, x in enumerate(Targetdomain[:, 1]) if x >= targetsite2])
    first = Targetdomain[indexa][0]
    last = Targetdomain[indexb][1]
    fragment = Targetori[first-1:last]
    
    # Create sequence file for UNAFold program
    target_seq_filename = f'userFiles/target_{access_code}.seq'
    with open(target_seq_filename, 'w') as file:
        file.write('>target\n')
        file.write(f'{fragment}\n')

    # Redefine fiveprime and threeprime
    fiveprime=targetsite1-first+1;
    threeprime=targetsite2-first+1;
    
    return fiveprime,threeprime


def deltaG3(fiveprime, threeprime, T, NA,access_code):
    # Build and execute command for mfold to calculate DG3f
    os.system(f"hybrid-ss-min userFiles/target_{access_code}.seq -n RNA -t {T} -T {T} -N {NA} -o userFiles/target_{access_code}")

    # Get the free energy from the mfold output
    DG3f = get_DG_from_file('userFiles/target_'+access_code+'.dG')

    # Calculate probe length
    length = threeprime - fiveprime + 1

    if DG3f>=max_folding_dG:
        DG3=max_folding_dG
    else:
        # Write constraints to the target.aux file
        with open('userFiles/target_'+access_code+'.aux', 'w') as file:
            file.write(f"P {fiveprime} 0 {length}\n")
            file.write(f"P 0 {fiveprime} {length}\n")

        # Build and execute command for mfold with constraints to calculate DG3u
        os.system(f"hybrid-ss-min userFiles/target_{access_code}.seq -n RNA -t {T} -T {T} -N {NA} -c --constraints=userFiles/target_{access_code}.aux -o userFiles/target_{access_code}")

        # Get the free energy from the mfold output with constraints
        DG3u = get_DG_from_file('userFiles/target_'+access_code+'.dG')

        if DG3u>=max_folding_dG:
            DG3=DG3f
        else:
            # Calculate final deltaG3
            DG3 = DG3f - DG3u

    return DG3, length

def deltaGov(T, DG1, DG2, DG3, Po=None):
    # Calculate RT
    Tk = 273.15 + T
    RT = R * Tk
    # Convert probe concentration from nM to M
    #Po = Po * 1E-9

    # Calculate Delta G overall
    DGov = -RT * math.log((math.exp(-DG1 / RT)) / ((math.exp(-DG2 / RT) + 1) * (math.exp(-DG3 / RT) + 1)))

    return DGov
    
def dGovHE(deltaGov, Po, T, access_code):
    # Constants and conversions
    Tk = 273.15 + T  # Convert Celsius to Kelvin
    RT = R * Tk
    Po = Po * 1E-9  # Convert probe concentration from nM to M

    # Generate x values for plotting
    x = np.arange(min(deltaGov-5, -20), 0.1, 0.1)
    y = deltaGov
    Kov = np.exp(-x / RT)
    hef = (Kov * Po) / (1 + Kov * Po)

    KovDG = np.exp(-y / RT)
    hefDG = (KovDG * Po) / (1 + KovDG * Po)

    # Define borders for plotting areas
    KovDGr = 1 / Po
    KovDGy = 999 / Po
    rightGREENborder = -RT * np.log(KovDGy)
    leftREDborder = -RT * np.log(KovDGr)
    GREENwidth = 4
    leftGREENborder = rightGREENborder - GREENwidth
    maxx = max(max(x), y) if deltaGov > max(x) else max(x)
    rightREDborder = maxx

    # Create the plot and save it
    f = plt.figure()
    Alpha=1.
    ax = f.add_subplot(111)
    ax.fill_between([leftGREENborder, rightGREENborder], 0, 1.2, color='green', alpha=Alpha)
    ax.fill_between([min(x), leftGREENborder], 0, 1.2, color='yellow', alpha=Alpha)
    ax.fill_between([rightGREENborder, leftREDborder], 0, 1.2, color=[249/256, 116/256, 28/256], alpha=Alpha)
    ax.fill_between([leftREDborder, maxx], 0, 1.2, color='red', alpha=Alpha)
    ax.plot(x, hef, 'k-', linewidth=2)
    ax.plot(y, hefDG, 'ms', markersize=10, label='This probe')
    ax.set_xlabel('ΔG (overall) [kcal/mol]')
    ax.set_ylabel('Hybridization Efficiency')
    ax.set_xlim([min(x), maxx])
    ax.set_ylim([0, 1.2])
    ax.legend()

    plt.savefig(f'userFiles/heplot_{access_code}.png')  # Save plot to file
    plt.close(f)

    return hefDG

def getHE(deltaGov, Po, T):
    # Constants and conversions
    Tk = 273.15 + T  # Convert Celsius to Kelvin
    RT = R * Tk
    Po = Po * 1E-9  # Convert probe concentration from nM to M

    Kov = np.exp(-deltaGov / RT)
    return (Kov * Po) / (1 + Kov * Po)


def fcurve(beta, x):
    teta1, teta2, teta5, m2, T, Po, L, dG1, dG2, dG3 = beta
    fa = x  # Formamide concentration percentage
    
    Po = Po * 1E-9 

    RT = 0.0019871 * (T + 273.15)  # Convert T from Celsius to Kelvin and calculate RT
    m1 = teta1 + teta2 * L
    m3 = teta5 * dG3

    K1 = np.exp(-(dG1 + m1 * fa) / RT)
    K2 = np.exp(-(dG2 + m2 * fa) / RT)
    K3 = np.exp(-(dG3 + m3 * fa) / RT)
    Kov = K1 / ((1 + K2) * (1 + K3))
    yhat = (Kov * Po) / (1 + Kov * Po)

    return yhat
    
def get_melting_point(x,y):
    if y[0]<0.5:
        return '< 0'
    elif min(y)>0.5:
        return '> '+str(int(max(x)))
    else:
        fam=np.nan
        i=1
        while i<len(x) and np.isnan(fam):
            if y[i]<0.5:
                dev=0.5-y[i]
                deltay=y[i-1]-y[i]
                deltax=x[i]-x[i-1]
                fam=round(x[i]-dev*deltax/deltay,1)
            i=i+1
        return fam

def formamideSmM(length, DG1, DG2, DG3, T, Po, access_code):
    teta1 = 0.0697
    teta2 = 0.0095
    teta5 = -0.0117
    m2 = 0.1

    beta = [teta1, teta2, teta5, m2, T, Po, length, DG1, DG2, DG3]
    fa = np.arange(0, 81)  # Formamide percentage from 0 to 80
    y = fcurve(beta, fa)
    
    #FAm
    fam=get_melting_point(fa,y)
    
    #Writing
    with open('userFiles/facurve_'+access_code+'.csv', 'w') as file:
        file.write('[FA] (%),Hyb. Eff.\n')
        for i in range(len(fa)):
            file.write(f"{fa[i]},{y[i]}\n")
    
    # Plotting
    plt.figure()
    plt.fill_between([0, max(fa)], 0, 1.2, color='yellow', alpha=0.5)  # Background
    plt.plot(fa, y, 'black', linewidth=2)  # Change color to match the attached figure
    plt.xlabel('Formamide Concentration (V/V)(%)')
    plt.ylabel('Hybridization Efficiency')
    plt.axis([min(fa), max(fa), 0, 1.2])
    plt.grid(True)
    plt.savefig('userFiles/faplot_'+access_code+'.png')  # Save the plot as a PNG file
    plt.close()

    return beta, fa, y, fam

def formamideSmM_dual(length, DG1, DG2, DG3, T, Po, access_code):
    teta1 = 0.0697
    teta2 = 0.0095
    teta5 = -0.0117
    m2 = 0.1

    # perfect match
    beta_pm = [teta1, teta2, teta5, m2, T, Po, length[0], DG1[0], DG2, DG3[0]]
    fa = np.arange(0, 81)  # Formamide percentage from 0 to 80
    y_pm = fcurve(beta_pm, fa)
    
    fam_pm=get_melting_point(fa,y_pm)
    
    # mismatch
    beta_mm = [teta1, teta2, teta5, m2, T, Po, length[1], DG1[1], DG2, DG3[1]]
    y_mm = fcurve(beta_mm, fa)
    
    fam_mm=get_melting_point(fa,y_mm)    
    
    
    #Writing
    with open('userFiles/facurve_'+access_code+'.csv', 'w') as file:
        file.write('[FA] (%),Hyb. Eff. (perfect match),Hyb. Eff. (mismatch)\n')
        for i in range(len(fa)):
            file.write(f"{fa[i]},{y_pm[i]},{y_mm[i]}\n")
    
    # Plotting
    plt.figure()
    plt.fill_between([0, max(fa)], 0, 1.2, color='yellow', alpha=0.5)  # Background
    plt.plot(fa, y_pm, 'black', linewidth=2, label='Target')
    plt.plot(fa, y_mm, 'green', linestyle='--', linewidth=2, label='Non-target')
    plt.xlabel('Formamide Concentration (V/V)(%)')
    plt.ylabel('Hybridization Efficiency')
    plt.axis([min(fa), max(fa), 0, 1.2])
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.savefig('userFiles/faplot_'+access_code+'.png')
    plt.close()

    return fam_pm,fam_mm
    
def fcurve_wcompetitor(beta_o, beta_1,beta_2, x):
    teta1, teta2, teta5, m2, T=beta_o
    Oo_1, L_1, dG1_1, dG2_1, dG3_1 = beta_1
    Oo_2, L_2, dG1_2, dG2_2, dG3_2 = beta_2
    fa = x  # Formamide concentration percentage
    
    Oo_1 = Oo_1 * 1E-9
    Oo_2 = Oo_2 * 1E-9

    RT = 0.0019871 * (T + 273.15)  # Convert T from Celsius to Kelvin and calculate RT
    
    m1_1 = teta1 + teta2 * L_1
    m1_2 = teta1 + teta2 * L_2
    
    m3_1 = teta5 * dG3_1
    m3_2 = teta5 * dG3_2

    K1_1 = np.exp(-(dG1_1 + m1_1 * fa) / RT)
    K1_2 = np.exp(-(dG1_2 + m1_2 * fa) / RT)
    
    K2_1 = np.exp(-(dG2_1 + m2 * fa) / RT)
    K2_2 = np.exp(-(dG2_2 + m2 * fa) / RT)
    
    K3_1 = np.exp(-(dG3_1 + m3_1 * fa) / RT)
    K3_2 = np.exp(-(dG3_2 + m3_2 * fa) / RT)
    
    Kov_1 = K1_1 / ((1 + K2_1) * (1 + K3_1))
    Kov_2 = K1_2 / ((1 + K2_2) * (1 + K3_2))
    
    yhat = (Kov_1 * Oo_1) / (1 + Kov_1 * Oo_1 + Kov_2 * Oo_2)

    return yhat
    
def formamideSmM_competitor(length, DG1, DG2, DG3, T, Oo, access_code):
    teta1 = 0.0697
    teta2 = 0.0095
    teta5 = -0.0117
    m2 = 0.1

    beta_o=[teta1, teta2, teta5, m2,T]
    fa = np.arange(0, 81)  # Formamide percentage from 0 to 80
    
    # Target
    beta_1 = [Oo[0], length[0], DG1[0], DG2[0], DG3[0]]
    beta_2 = [Oo[2], length[2], DG1[2], DG2[2], DG3[2]]    
    
    y_target = fcurve_wcompetitor(beta_o,beta_1,beta_2,fa)
    
    # Nontarget
    beta_1 = [Oo[1], length[1], DG1[1], DG2[1], DG3[1]]
    beta_2 = [Oo[3], length[3], DG1[3], DG2[3], DG3[3]]    
    
    y_nontarget = fcurve_wcompetitor(beta_o,beta_1,beta_2,fa)   
    
    #Writing
    with open('userFiles/facurve_'+access_code+'.csv', 'w') as file:
        file.write('[FA] (%),Hyb. Eff. (perfect match),Hyb. Eff. (mismatch)\n')
        for i in range(len(fa)):
            file.write(f"{fa[i]},{y_target[i]},{y_nontarget[i]}\n")
    
    # Plotting
    plt.figure()
    plt.fill_between([0, max(fa)], 0, 1.2, color='yellow', alpha=0.5)  # Background
    plt.plot(fa, y_target, 'black', linewidth=2, label='Target')
    plt.plot(fa, y_nontarget, 'green', linestyle='--', linewidth=2, label='Non-target')
    plt.xlabel('Formamide Concentration (V/V)(%)')
    plt.ylabel('Hybridization Efficiency')
    plt.axis([min(fa), max(fa), 0, 1.2])
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.savefig('userFiles/faplot_'+access_code+'.png')
    plt.close()


def parse_dnd(accode,L,name):
    '''Parses out the first distance value in dnd files from ClustalW alignments. If greater than the tolerable value, adds to list <L> as a tuple with <name>'''
    with open('userFiles/combined_'+accode+'.dnd') as f:
        line = f.readline()
    match = re.search(r'R:([0-9.]+)', line)
    dist=float(match.group(1)) if match else np.nan
    if np.isnan(dist) or dist>maxDist:
        L.append((name,dist))

def get_alignment_warning(bad_alignments,targetmolecule,domain,help_anchor,higher_folding_dG=-9999):
    '''Returns either a warning for large distances observed in aligned pairs or a reference to a help section.'''
    if higher_folding_dG>=max_folding_dG:
        s_folding=f'<br><br> <b>Note</b>: Folding free energies are capped at {max_folding_dG} kcal/mol. See help documentation for free energy changes of the <a href="helpdocuments.html#deltaG2" target="_blank">second</a> and <a href="helpdocuments.html#deltaG3" target="_blank">third</a> reactions.'
    else:
        s_folding=''
    if bad_alignments:
        if len(bad_alignments)>1:
            sfx,a='s',''
        else:
            sfx,a='','a '
        s1='<b><u>WARNING: Bad seqeuence alignment! </b></u><br>'
        s2='During the calculation of deltaG3, the targeted sequences are aligned with reference sequences. '
        s3=f'The following alignment{sfx} resulted in {a}large distance{sfx}:<br><br>'
        s4='<br>'.join([f'{t[0]} sequence against the {targetmolecule} rRNA of {ref_organisms[domain][0]} reference <i>{ref_organisms[domain][1]}</i> (distance={t[1]:.2f})' for t in bad_alignments])
        s5=f'<br><br> Please make sure the relevant entries (target molecule=<b>{targetmolecule}</b>, domain=<b>{ref_organisms[domain][0]}</b>) are correct and redo the calculation if necessary.'
        s6=f'<br><br> See <a href="helpdocuments.html#{help_anchor}" target="_blank">help</a> for more information.'+s_folding
        return s1+s2+s3+s4+s5+s6
    else:
        return f'See <a href="helpdocuments.html#{help_anchor}" target="_blank">help</a> for the interpretation of results.'+s_folding
        
def formamideSmM_multiple(Length, DG1, DG2, DG3, T, Po, access_code):
    teta1 = 0.0697
    teta2 = 0.0095
    teta5 = -0.0117
    m2 = 0.1
    
    fa = np.arange(0, 81)  # Formamide percentage from 0 to 80
    Y=[]
    for i in range(len(DG1)):
        beta = [teta1, teta2, teta5, m2, T, Po, Length[i], DG1[i], DG2[i], DG3[i]]
        y = fcurve(beta, fa)
        Y.append(y)    
    
    #Writing
    with open('userFiles/facurve_'+access_code+'.csv', 'w') as file:
        file.write('[FA] (%),'+','.join(['Hyp. Eff. (Target-'+str(i+1)+')' for i in range(len(Length))])+'\n')
        for i in range(len(fa)):
            file.write(f"{fa[i]},{','.join([str(y[i]) for y in Y])}\n")
    
    # Plotting
    colors=['black','red','green','orange','magenta']
    plt.figure()
    plt.fill_between([0, max(fa)], 0, 1.2, color='yellow', alpha=0.5)  # Background
    for i in range(len(Y)):
        plt.plot(fa, Y[i], colors[i], linewidth=2, label=str(i+1))
    plt.xlabel('Formamide Concentration (V/V)(%)')
    plt.ylabel('Hybridization Efficiency')
    plt.axis([min(fa), max(fa), 0, 1.2])
    plt.grid(True)
    plt.legend(loc='upper left',ncol=len(Y))
    plt.savefig('userFiles/faplot_'+access_code+'.png')
    plt.close()      
        
        