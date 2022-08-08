#version 2.2 - 26/05/2020

# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import time

"""It combines the pickle binary data produced by the get_info.py script."""


#we get the name of input and output files
N_argument_files=len(sys.argv)-1

if(N_argument_files<3):
   print('Syntax: ./combine.py <output pickle binary data file> <input pickle binary data 1> [input pickle binary data 2] ...')
   print('(Please, note that suffixes will be automatically added to the two output files, but the names of the input files must be complete.)')
   sys.exit(1)

outputfile=sys.argv[1]
inputfiles=sys.argv[2:]
nif=len(inputfiles)

thist=[]
xhist=[]
yhist=[]
zhist=[]
thist_cp=[]
xhist_cp=[]
yhist_cp=[]
zhist_cp=[]
temphist=[]
muhist=[]
betahist=[]
phase=[]
Tmupoints=[]
Tpt=[]
mupt=[]
muoverTpt=[]
Trap=[]
murap=[]
muoverTrap=[]
pcount=[]
pcount_raplim=[]
ntothist=[]
e_over_nhist=[]
sT3hist=[]

first_time=True
for i in range(0,nif):
  if(not os.path.isfile(inputfiles[i])):
    print("Warning, input file "+inputfiles[i]+" does not exist!")
  else:
    print("Reading file "+inputfiles[i])
    with open(inputfiles[i],"rb") as pi:
          data=pickle.load(pi)
    if(first_time):
      ptypeF,pcountF,pcount_raplimF,raplimF,ptlimF,mbF,tbF,ptarrF,raparrF,thist_bigF,xhist_bigF,yhist_bigF,zhist_bigF,thist_big_ovF,xhist_big_ovF,yhist_big_ovF,zhist_big_ovF,thist_cpF,xhist_cpF,yhist_cpF,zhist_cpF,thistF,xhistF,yhistF,zhistF,temphistF,muhistF,betahistF,phaseF,TmupointsF,TptF,muptF,muoverTptF,TrapF,murapF,muoverTrapF,ntothistF,e_over_nhistF,sT3histF=data[:]       
      nitem=len(thistF)
      for j in range(0,nitem):
        thist.append(np.asarray(thistF[j][0]))
        xhist.append(np.asarray(xhistF[j][0]))
        yhist.append(np.asarray(yhistF[j][0]))
        zhist.append(np.asarray(zhistF[j][0]))
        thist_cp.append(np.asarray(thist_cpF[j][0]))
        xhist_cp.append(np.asarray(xhist_cpF[j][0]))
        yhist_cp.append(np.asarray(yhist_cpF[j][0]))
        zhist_cp.append(np.asarray(zhist_cpF[j][0]))
        temphist.append(np.asarray(temphistF[j][0]))
        muhist.append(np.asarray(muhistF[j][0]))
        betahist.append(np.asarray(betahistF[j][0]))
        ntothist.append(np.asarray(ntothistF[j][0]))
        e_over_nhist.append(np.asarray(e_over_nhistF[j][0]))
        sT3hist.append(np.asarray(sT3histF[j][0]))
        phase.append(np.asarray(phaseF[j][0]))
        Tmupoints.append(np.asarray(TmupointsF[j][0:]))
        Tpt.append(np.asarray(TptF[j][0:]))
        mupt.append(np.asarray(muptF[j][0:]))
        muoverTpt.append(np.asarray(muoverTptF[j][0:]))
        Trap.append(np.asarray(TrapF[j][0:]))
        murap.append(np.asarray(murapF[j][0:]))
        muoverTrap.append(np.asarray(muoverTrapF[j][0:]))
      pcount=np.asarray(pcountF)
      pcount_raplim=np.asarray(pcount_raplimF)
      thist_big=np.asarray(thist_bigF[0])
      xhist_big=np.asarray(xhist_bigF[0])
      yhist_big=np.asarray(yhist_bigF[0])
      zhist_big=np.asarray(zhist_bigF[0])
      thist_big_ov=np.asarray(thist_big_ovF[0])
      xhist_big_ov=np.asarray(xhist_big_ovF[0])
      yhist_big_ov=np.asarray(yhist_big_ovF[0])
      zhist_big_ov=np.asarray(zhist_big_ovF[0])
      raplim=raplimF
      ptlim=ptlimF
      first_time=False
    else:
      ptypeN,pcountN,pcount_raplimN,raplimN,ptlimN,mbN,tbN,ptarrN,raparrN,thist_bigN,xhist_bigN,yhist_bigN,zhist_bigN,thist_big_ovN,xhist_big_ovN,yhist_big_ovN,zhist_big_ovN,thist_cpN,xhist_cpN,yhist_cpN,zhist_cpN,thistN,xhistN,yhistN,zhistN,temphistN,muhistN,betahistN,phaseN,TmupointsN,TptN,muptN,muoverTptN,TrapN,murapN,muoverTrapN,ntothistN,e_over_nhistN,sT3histN=data[:]      
      if(len(thistN)!=nitem):
        print("Error, it seems that the number of particle types has changed, because the new thist array of histograms has a different number of elements..")
        print("Are you sure that you are combining the files correctly?")
        sys.exit(2)
      pcount+=np.asarray(pcountN)
      pcount_raplim+=np.asarray(pcount_raplimN)
      if(ptlimN != ptlim):
        print("Error, there is a mismatching between the pT range of the previous file(s) and the current one: "+str(ptlim)+" vs "+str(ptlimN)+".")
        print("Sorry, but combining them might lead to inconsistencies. I prefer to quit.")
        sys.exit(2)
      if(raplimN != raplim):
        print("Error, there is a mismatching between the rapidity range of the previous file(s) and the current one: "+str(raplim)+" vs "+str(raplimN)+".")
        print("Sorry, but combining them might lead to inconsistencies. I prefer to quit.")
        sys.exit(2)
      for j in range(0,nitem):
        new_data=np.asarray(thistN[j][0])
        thist[j]+=new_data
        new_data=np.asarray(xhistN[j][0])
        xhist[j]+=new_data
        new_data=np.asarray(yhistN[j][0])
        yhist[j]+=new_data
        new_data=np.asarray(zhistN[j][0])
        zhist[j]+=new_data

        new_data=np.asarray(thist_cpN[j][0])
        thist_cp[j]+=new_data
        new_data=np.asarray(xhist_cpN[j][0])
        xhist_cp[j]+=new_data
        new_data=np.asarray(yhist_cpN[j][0])
        yhist_cp[j]+=new_data
        new_data=np.asarray(zhist_cpN[j][0])
        zhist_cp[j]+=new_data

        new_data=np.asarray(temphistN[j][0])
        temphist[j]+=new_data
        new_data=np.asarray(muhistN[j][0])
        muhist[j]+=new_data
        new_data=np.asarray(betahistN[j][0])
        betahist[j]+=new_data
        new_data=np.asarray(ntothistN[j][0])
        ntothist[j]+=new_data
        new_data=np.asarray(e_over_nhistN[j][0])
        e_over_nhist[j]+=new_data
        new_data=np.asarray(sT3histN[j][0])
        sT3hist[j]+=new_data
        new_data=np.asarray(phaseN[j][0])
        phase[j]+=new_data
        TmupC,Tmup_errC,muC,mu_errC,betaC,beta_errC,pointnumC=TmupointsN[j][0:]
        nold=Tmupoints[j][-1]
        T_avgold=Tmupoints[j][0]
        mu_avgold=Tmupoints[j][2]
        beta_avgold=Tmupoints[j][4]
        pointnumR=nold+pointnumC
        if(pointnumR != 0):
          Tmupoints[j][0]=(nold*T_avgold+TmupC*pointnumC)/pointnumR
          delta1=Tmupoints[j][0]-T_avgold
          delta2=Tmupoints[j][0]-TmupC
        else:
          Tmupoints[j][0]=0.
        if(pointnumR > 1):
          Tmupoints[j][1]=math.sqrt(((Tmupoints[j][1]**2)*(nold-1)+(Tmup_errC**2)*(pointnumC-1)+nold*delta1**2+pointnumC*delta2**2)/(pointnumR-1))
        else:
          Tmupoints[j][1]=0.
        if(pointnumR != 0):
          Tmupoints[j][2]=(nold*mu_avgold+muC*pointnumC)/pointnumR
          delta1=Tmupoints[j][2]-mu_avgold
          delta2=Tmupoints[j][2]-muC
        else:
          Tmupoints[j][2]=0.
        if(pointnumR > 1):
          Tmupoints[j][3]=math.sqrt(((Tmupoints[j][3]**2)*(nold-1)+(mu_errC**2)*(pointnumC-1)+nold*delta1**2+pointnumC*delta2**2)/(pointnumR-1))
        else:
          Tmupoints[j][3]=0.
        if(pointnumR != 0):
          Tmupoints[j][4]=(nold*beta_avgold+betaC*pointnumC)/pointnumR
          delta1=Tmupoints[j][4]-beta_avgold
          delta2=Tmupoints[j][4]-betaC
        else:
          Tmupoints[j][4]=0.
        if(pointnumR > 1):
          Tmupoints[j][5]=math.sqrt(((Tmupoints[j][5]**2)*(nold-1)+(beta_errC**2)*(pointnumC-1)+nold*delta1**2+pointnumC*delta2**2)/(pointnumR-1))
        else:
          Tmupoints[j][5]=0.
        Tmupoints[j][6]=pointnumR
        
        Tpt_avg,Tpt_err,Tpt_count=TptN[j][0:]
        Tpt_count_old=Tpt[j][-1]
        Tpt_newcount=Tpt_count+Tpt_count_old
        for k in range(len(Tpt_newcount)):
          if(Tpt_newcount[k] != 0):
            avgold=Tpt[j][0][k]
            Tpt[j][0][k]=(Tpt_count_old[k]*avgold+Tpt_avg[k]*Tpt_count[k])/Tpt_newcount[k]
            delta1=Tpt[j][0][k]-avgold
            delta2=Tpt[j][0][k]-Tpt_avg[k]
          else:
            Tpt[j][0][k]=0.
          if(Tpt_newcount[k] > 1):
            Tpt[j][1][k]=math.sqrt(((Tpt[j][1][k]**2)*(Tpt_count_old[k]-1)+(Tpt_err[k]**2)*(Tpt_count[k]-1)+Tpt_count_old[k]*delta1**2+Tpt_count[k]*delta2**2)/(Tpt_newcount[k]-1))
          else:
            Tpt[j][1][k]=0.
          Tpt[j][2][k]=Tpt_newcount[k]
       
        mupt_avg,mupt_err,mupt_count=muptN[j][0:]
        mupt_count_old=mupt[j][-1]
        mupt_newcount=mupt_count+mupt_count_old
        for k in range(len(mupt_newcount)):
          if(mupt_newcount[k] != 0):
            avgold=mupt[j][0][k]
            mupt[j][0][k]=(mupt_count_old[k]*avgold+mupt_avg[k]*mupt_count[k])/mupt_newcount[k]
            delta1=mupt[j][0][k]-avgold
            delta2=mupt[j][0][k]-mupt_avg[k]
          else:
            mupt[j][0][k]=0.
          if(mupt_newcount[k] > 1):
            mupt[j][1][k]=math.sqrt(((mupt[j][1][k]**2)*(mupt_count_old[k]-1)+(mupt_err[k]**2)*(mupt_count[k]-1)+mupt_count_old[k]*delta1**2+mupt_count[k]*delta2**2)/(mupt_newcount[k]-1))
          else:
            mupt[j][1][k]=0.
          mupt[j][2][k]=mupt_newcount[k]

        muoverTpt_avg,muoverTpt_err,muoverTpt_count=muoverTptN[j][0:]
        muoverTpt_count_old=muoverTpt[j][-1]
        muoverTpt_newcount=muoverTpt_count+muoverTpt_count_old
        for k in range(len(muoverTpt_newcount)):
          if(muoverTpt_newcount[k] != 0):
            avgold=muoverTpt[j][0][k]
            muoverTpt[j][0][k]=(muoverTpt_count_old[k]*avgold+muoverTpt_avg[k]*muoverTpt_count[k])/muoverTpt_newcount[k]
            delta1=muoverTpt[j][0][k]-avgold
            delta2=muoverTpt[j][0][k]-muoverTpt_avg[k]
          else:
            muoverTpt[j][0][k]=0.
          if(muoverTpt_newcount[k] > 1):
            muoverTpt[j][1][k]=math.sqrt(((muoverTpt[j][1][k]**2)*(muoverTpt_count_old[k]-1)+(muoverTpt_err[k]**2)*(muoverTpt_count[k]-1)+muoverTpt_count_old[k]*delta1**2+muoverTpt_count[k]*delta2**2)/(muoverTpt_newcount[k]-1))
          else:
            muoverTpt[j][1][k]=0.
          muoverTpt[j][2][k]=muoverTpt_newcount[k]

        murap_avg,murap_err,murap_count=murapN[j][0:]
        murap_count_old=murap[j][-1]
        murap_newcount=murap_count+murap_count_old
        for k in range(len(murap_newcount)):
          if(murap_newcount[k] != 0):
            avgold=murap[j][0][k]
            murap[j][0][k]=(murap_count_old[k]*avgold+murap_avg[k]*murap_count[k])/murap_newcount[k]
            delta1=murap[j][0][k]-avgold
            delta2=murap[j][0][k]-murap_avg[k]
          else:
            murap[j][0][k]=0.
          if(murap_newcount[k] > 1):
            murap[j][1][k]=math.sqrt(((murap[j][1][k]**2)*(murap_count_old[k]-1)+(murap_err[k]**2)*(murap_count[k]-1)+murap_count_old[k]*delta1**2+murap_count[k]*delta2**2)/(murap_newcount[k]-1))
          else:
            murap[j][1][k]=0.
          murap[j][2][k]=murap_newcount[k]

        muoverTrap_avg,muoverTrap_err,muoverTrap_count=muoverTrapN[j][0:]
        muoverTrap_count_old=muoverTrap[j][-1]
        muoverTrap_newcount=muoverTrap_count+muoverTrap_count_old
        for k in range(len(muoverTrap_newcount)):
          if(muoverTrap_newcount[k] != 0):
            avgold=muoverTrap[j][0][k]
            muoverTrap[j][0][k]=(muoverTrap_count_old[k]*avgold+muoverTrap_avg[k]*muoverTrap_count[k])/muoverTrap_newcount[k]
            delta1=muoverTrap[j][0][k]-avgold
            delta2=muoverTrap[j][0][k]-muoverTrap_avg[k]
          else:
            muoverTrap[j][0][k]=0.
          if(muoverTrap_newcount[k] > 1):
            muoverTrap[j][1][k]=math.sqrt(((muoverTrap[j][1][k]**2)*(muoverTrap_count_old[k]-1)+(muoverTrap_err[k]**2)*(muoverTrap_count[k]-1)+muoverTrap_count_old[k]*delta1**2+muoverTrap_count[k]*delta2**2)/(muoverTrap_newcount[k]-1))
          else:
            muoverTrap[j][1][k]=0.
          muoverTrap[j][2][k]=muoverTrap_newcount[k]

        Trap_avg,Trap_err,Trap_count=TrapN[j][0:]
        Trap_count_old=Trap[j][-1]
        Trap_newcount=Trap_count+Trap_count_old
        for k in range(len(Trap_newcount)):
          if(Trap_newcount[k] != 0):
            avgold=Trap[j][0][k]
            Trap[j][0][k]=(Trap_count_old[k]*avgold+Trap_avg[k]*Trap_count[k])/Trap_newcount[k]
            delta1=Trap[j][0][k]-avgold
            delta2=Trap[j][0][k]-Trap_avg[k]
          else:
            Trap[j][0][k]=0.
          if(Trap_newcount[k] > 1):
            Trap[j][1][k]=math.sqrt(((Trap[j][1][k]**2)*(Trap_count_old[k]-1)+(Trap_err[k]**2)*(Trap_count[k]-1)+Trap_count_old[k]*delta1**2+Trap_count[k]*delta2**2)/(Trap_newcount[k]-1))
          else:
            Trap[j][1][k]=0.
          Trap[j][2][k]=Trap_newcount[k]

      new_data=np.asarray(thist_bigN[0])
      thist_big+=new_data
      new_data=np.asarray(xhist_bigN[0])
      xhist_big+=new_data
      new_data=np.asarray(yhist_bigN[0])
      yhist_big+=new_data
      new_data=np.asarray(zhist_bigN[0])
      zhist_big+=new_data
      new_data=np.asarray(thist_big_ovN[0])
      thist_big_ov+=new_data
      new_data=np.asarray(xhist_big_ovN[0])
      xhist_big_ov+=new_data
      new_data=np.asarray(yhist_big_ovN[0])
      yhist_big_ov+=new_data
      new_data=np.asarray(zhist_big_ovN[0])
      zhist_big_ov+=new_data

print("Reassembling histograms as tuples")
ptypeE=ptypeF
pcountE=pcount
pcount_raplimE=pcount_raplim
mbE=mbF
tbE=tbF
ptarrE=ptarrF
raparrE=raparrF
thist_bigE=(thist_big,thist_bigF[1])
xhist_bigE=(xhist_big,xhist_bigF[1])
yhist_bigE=(yhist_big,yhist_bigF[1])
zhist_bigE=(zhist_big,zhist_bigF[1])
thist_big_ovE=(thist_big_ov,thist_big_ovF[1])
xhist_big_ovE=(xhist_big_ov,xhist_big_ovF[1])
yhist_big_ovE=(yhist_big_ov,yhist_big_ovF[1])
zhist_big_ovE=(zhist_big_ov,zhist_big_ovF[1])
thistE=[]
xhistE=[]
yhistE=[]
zhistE=[]
thist_cpE=[]
xhist_cpE=[]
yhist_cpE=[]
zhist_cpE=[]
temphistE=[]
muhistE=[]
betahistE=[]
ntothistE=[]
e_over_nhistE=[]
sT3histE=[]
phaseE=[]
muptE=[]
muoverTptE=[]
TptE=[]
murapE=[]
muoverTrapE=[]
TrapE=[]
TmupointsE=[]
for i in range(0,nitem):
    thistE.append((thist[i],thistF[i][1]))
    xhistE.append((xhist[i],xhistF[i][1]))
    yhistE.append((yhist[i],yhistF[i][1]))
    zhistE.append((zhist[i],zhistF[i][1]))
    thist_cpE.append((thist_cp[i],thistF[i][1]))
    xhist_cpE.append((xhist_cp[i],xhistF[i][1]))
    yhist_cpE.append((yhist_cp[i],yhistF[i][1]))
    zhist_cpE.append((zhist_cp[i],zhistF[i][1]))
    temphistE.append((temphist[i],temphistF[i][1]))
    muhistE.append((muhist[i],muhistF[i][1]))
    betahistE.append((betahist[i],betahistF[i][1]))
    ntothistE.append((ntothist[i],ntothistF[i][1]))
    e_over_nhistE.append((e_over_nhist[i],e_over_nhistF[i][1]))
    sT3histE.append((sT3hist[i],sT3histF[i][1]))
    phaseE.append((phase[i],phaseF[i][1],phaseF[i][2]))
    TptE.append((Tpt[i][0],Tpt[i][1],Tpt[i][2]))
    muptE.append((mupt[i][0],mupt[i][1],mupt[i][2]))
    muoverTptE.append((muoverTpt[i][0],muoverTpt[i][1],muoverTpt[i][2]))
    TrapE.append((Trap[i][0],Trap[i][1],Trap[i][2]))
    murapE.append((murap[i][0],murap[i][1],murap[i][2]))
    muoverTrapE.append((muoverTrap[i][0],muoverTrap[i][1],muoverTrap[i][2]))
    TmupointsE.append((Tmupoints[i][0],Tmupoints[i][1],Tmupoints[i][2],Tmupoints[i][3],Tmupoints[i][4],Tmupoints[i][5],Tmupoints[i][6]))

print("Pickling ptype,pcount,pcount_raplim,raplim,ptlim,mb,tb,ptarr,raparr,thist_big,xhist_big,yhist_big,zhist_big,thist_big_ov,xhist_big_ov,yhist_big_ov,zhist_big_ov,thist,xhist,yhist,zhist,temphist,muhist,phase,Tmupoints,Tpt,mupt,Trap,murap as tuples of histograms")
pickout=outputfile+"_plot_data_binary.pickle"
if(os.path.isfile(pickout)):
  pickout_old=pickout
  timestring=str(int(time.time()))
  pickout=pickout+"_"+timestring
  print("Output file "+pickout_old+" already exists, renaming pickle output file as: "+pickout)
with open(pickout,"wb") as po:
     pickle.dump((ptypeE,pcountE,pcount_raplimE,raplim,ptlim,mbE,tbE,ptarrE,raparrE,thist_bigE,xhist_bigE,yhist_bigE,zhist_bigE,thist_big_ovE,xhist_big_ovE,yhist_big_ovE,zhist_big_ovE,thist_cpE,xhist_cpE,yhist_cpE,zhist_cpE,thistE,xhistE,yhistE,zhistE,temphistE,muhistE,betahistE,phaseE,TmupointsE,TptE,muptE,muoverTptE,TrapE,murapE,muoverTrapE,ntothistE,e_over_nhistE,sT3histE),po)
print("Printing text data for plots with results from different simulations")
sp="    "
dout=outputfile+"_point_data.txt"
if(os.path.isfile(dout)):
  dout_old=dout
  dout=dout+"_"+timestring
  print("Output file "+dout_old+" already exists, renaming pickle output file as: "+dout)
po=open(dout,"w")
po.write("#Particle    T   T_error    mu     mu_error     <beta>     <beta>_err\n")
for i,v in ptypeE.items():
    indx=v[0]
    po.write(v[1]+sp+'{:7.3f}'.format(TmupointsE[indx][0])+sp+'{:7.3f}'.format(TmupointsE[indx][1])+sp+'{:7.3f}'.format(TmupointsE[indx][2])+sp+'{:7.3f}'.format(TmupointsE[indx][3])+sp+'{:7.3f}'.format(TmupointsE[indx][4])+sp+'{:7.3f}'.format(TmupointsE[indx][5])+"\n")
po.write("All"+sp+'{:7.3f}'.format(TmupointsE[-1][0])+sp+'{:7.3f}'.format(TmupointsE[-1][1])+sp+'{:7.3f}'.format(TmupointsE[-1][2])+sp+'{:7.3f}'.format(TmupointsE[-1][3])+sp+'{:7.3f}'.format(TmupointsE[-1][4])+sp+'{:7.3f}'.format(TmupointsE[-1][5])+"\n")
po.close()
print("All done.")
