from utils import *

def dG1_series(vars,accode,Targetori,ts,L):
    ## Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])
    stack, deltaH, deltaS = read_stacking_file('inputFiles/stacking.txt')    
    
    ## dG1
    DG1,Probe,TS=[],[],[]
    ts1=ts[0]
    ts2=ts1+L-1
    while ts2<=ts[1]:
        probe=fprobeseq(Targetori,ts1,ts2);
        Probe.append(probe)
        TS.append((ts1,ts2))
        DG1.append(fnearestneighbor(probe,stack,deltaH,deltaS,T,NA))
        ts1=ts1+1
        ts2=ts1+L-1


    return DG1, TS, Probe
    
def dG2_series(vars,accode,Targetori,ts,L):
    ## Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])

    ## dG2
    DG2,Probe,TS=[],[],[]
    ts1=ts[0]
    ts2=ts1+L-1
    while ts2<=ts[1]:
        probe=fprobeseq(Targetori,ts1,ts2);
        Probe.append(probe)
        TS.append((ts1,ts2))
        DG2.append(deltaG2(probe, T, NA, accode))
        ts1=ts1+1
        ts2=ts1+L-1
    
    return DG2, TS, Probe
    
def dG3_series(vars,accode,Targetori,ts,L):
    ## Thermo
    T=float(vars['T'])
    NA=float(vars['NA'])
    targetmolecule=vars['targetmolecule']
    domain=vars['domain']
    
    #dG3
    ts1=ts[0]
    ts2=ts1+L-1   
    if len(targetmolecule)==3:
        fiveprime,threeprime=newaligner(Targetori, ts1, ts2, targetmolecule, domain, accode);
    else:
        with open('userFiles/target_'+accode+'.seq', 'w') as file:
            file.write('>target\n')
            file.write(f'{Targetori}\n')
            fiveprime,threeprime=(ts1,ts2)
    #print(targetsite,fiveprime,threeprime)
    DG3_first, _ = deltaG3(fiveprime, threeprime, T, NA, accode);
    #print(fiveprime,threeprime,DG3_first)
    
    DG3,Probe,TS=[DG3_first],[fprobeseq(Targetori,ts1,ts2)],[(ts1,ts2)]
    ts1=ts1+1
    ts2=ts1+L-1
    while ts2<=ts[1]:
        probe=fprobeseq(Targetori,ts1,ts2);
        Probe.append(probe)
        TS.append((ts1,ts2))
        if len(targetmolecule)==3:
            fiveprime,threeprime=newaligner(Targetori, ts1, ts2, targetmolecule, domain, accode,skip_alignment=True);
        else:
            fiveprime,threeprime=(ts1,ts2)
        DG3_this,_=deltaG3(fiveprime, threeprime, T, NA, accode);
        DG3.append(DG3_this)
        ts1=ts1+1
        ts2=ts1+L-1
    
    return DG3, TS, Probe

def writeAndPlot(DG,TS,Probe,accode,n):
        #Writing
        with open('userFiles/dg'+str(n)+'series_'+accode+'.csv', 'w') as file:
            file.write("5' end position,3' end position,Probe sequence,deltaG"+str(n)+" (kcal/mol)\n")
            for i in range(len(DG)):
                file.write(f"{TS[i][0]},{TS[i][1]},{Probe[i]},{DG[i]}\n")
        
        # Plotting
        x = [str(int(ts[0])) for ts in TS]
        y = [round(val,1) for val in DG]
        plt.figure()
        
        plt.plot(x, y, color='orange', linewidth=1)
        plt.scatter(x, y, color='red', s=20)
        
        plt.xlabel("5' end", fontweight='bold')
        plt.ylabel("ΔG°"+str(n)+" (kcal/mol)", fontweight='bold')
        
        step = 1 if len(x) <= 25 else 2  
        tick_positions = range(0, len(x), step)
        tick_labels = [x[i] for i in tick_positions]
        plt.xticks(ticks=tick_positions, labels=tick_labels, rotation=90)

        plt.grid(True, axis='both')  # ensure grid matches all ticks
        plt.tight_layout()
        plt.savefig('userFiles/dg'+str(n)+'plot_'+accode+'.png')
        plt.close()        
    
    