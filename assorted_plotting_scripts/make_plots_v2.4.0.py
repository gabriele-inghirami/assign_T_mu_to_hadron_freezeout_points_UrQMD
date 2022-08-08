# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3
#version 2.4.0 - 26/5/2020
#it goes with get_info_v2.2.0.py written in May 2020

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import matplotlib
import matplotlib.pyplot as plt

"""It creates plots from the pickle binary data produced by the get_info.py script."""


#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files!=2):
   print ('Syntax: ./make_plots.py <pickle binary data file> <output dir>')
   sys.exit(1)

inputfile=sys.argv[1]

od=sys.argv[2]
if(not os.path.exists(od)):
  os.mkdir(od)


with open(inputfile,"rb") as pi:
     data=pickle.load(pi)

ptype,pcount,pcount_raplim,raplim,ptlim,mb,tb,ptarr,raparr,thist_big,xhist_big,yhist_big,zhist_big,thist_big_ov,xhist_big_ov,yhist_big_ov,zhist_big_ov,thist_cp,xhist_cp,yhist_cp,zhist_cp,thist,xhist,yhist,zhist,temphist,muhist,betahist,phase,Tmupoints,Tpt,mupt,muoverTpt,Trap,murap,muoverTrap,ntot,e_over_n,sT3=data[:]

symb={"p":r"$p$","N":r"$n$","ap":r"$\bar{p}$","aN":r"$\bar{n}$","L":r"$\Lambda","aL":r"$\bar \Lambda$","S0":r"\Sigma^0$","aS0":r"$\bar{\Sigma}^0$","Xi0":r"$\Xi^0$","aXi0":r"$\bar{\Xi}^0$","Om":r"$\Omega$","aOm":r"$\bar \Omega$","pip":r"$\pi^+$","pi0":r"$\pi^0$","pim":r"$\pi^-$","Kp":r"$K^+$","Km":r"$K^-$","K0":r"$K^0$","aK0":r"$\bar{K}^0$","deu":"d","adeu":r"$\bar{d}$"}

pnames=[""]*len(pcount)
psymb=[""]*len(pcount)
for i,v in ptype.items():
    pnames[v[0]]=v[1]
    psymb[v[0]]=symb[v[1]]
pnames[-1]="All"  
psymb[-1]="All"  

#we assume that all histograms and plots have the same axes, regardless of the particle specie
lt=len(thist[-1][1])-1
tt=np.zeros(lt)
lx=len(xhist[-1][1])-1
xx=np.zeros(lx)
ly=len(yhist[-1][1])-1
yy=np.zeros(ly)
lz=len(zhist[-1][1])-1
zz=np.zeros(lz)
lte=len(temphist[-1][1])-1
ate=np.zeros(lte)
lm=len(muhist[-1][1])-1
am=np.zeros(lm)
lntot=len(ntot[-1][1])-1
antot=np.zeros(lntot)
leovern=len(e_over_n[-1][1])-1
aeovern=np.zeros(leovern)
lsT3=len(sT3[-1][1])-1
asT3=np.zeros(lsT3)
lb=len(betahist[-1][1])-1
ab=np.zeros(lb)
lpt=len(ptarr)-1
apt=np.zeros(lpt)
lrap=len(raparr)-1
arap=np.zeros(lrap)
for i in range(lt):
  tt[i]=(thist[-1][1][i+1]+thist[-1][1][i])/2.
for i in range(lx):
  xx[i]=(xhist[-1][1][i+1]+xhist[-1][1][i])/2.
for i in range(ly):
  yy[i]=(yhist[-1][1][i+1]+yhist[-1][1][i])/2.
for i in range(lz):
  zz[i]=(zhist[-1][1][i+1]+zhist[-1][1][i])/2.
for i in range(lm):
  am[i]=(muhist[-1][1][i+1]+muhist[-1][1][i])/2.
for i in range(lb):
  ab[i]=(betahist[-1][1][i+1]+betahist[-1][1][i])/2.
for i in range(lte):
  ate[i]=(temphist[-1][1][i+1]+temphist[-1][1][i])/2.
for i in range(lntot):
  antot[i]=(ntot[-1][1][i+1]+ntot[-1][1][i])/2.
for i in range(leovern):
  aeovern[i]=(e_over_n[-1][1][i+1]+e_over_n[-1][1][i])/2.
for i in range(lsT3):
  asT3[i]=(sT3[-1][1][i+1]+sT3[-1][1][i])/2.
for i in range(lpt):
  apt[i]=(ptarr[i+1]+ptarr[i])/2.
for i in range(lrap):
  arap[i]=(raparr[i+1]+raparr[i])/2.
#same as before, but now for the whole range of particles before associating coarse graining data
lt_ov=len(thist_big_ov[1])-1
tt_ov=np.zeros(lt_ov)
lx_ov=len(xhist_big_ov[1])-1
xx_ov=np.zeros(lx_ov)
ly_ov=len(yhist_big_ov[1])-1
yy_ov=np.zeros(ly_ov)
lz_ov=len(zhist_big_ov[1])-1
zz_ov=np.zeros(lz_ov)
for i in range(lt_ov):
  tt_ov[i]=(thist_big_ov[1][i+1]+thist_big_ov[1][i])/2.
for i in range(lx_ov):
  xx_ov[i]=(xhist_big_ov[1][i+1]+xhist_big_ov[1][i])/2.
for i in range(ly_ov):
  yy_ov[i]=(yhist_big_ov[1][i+1]+yhist_big_ov[1][i])/2.
for i in range(lz_ov):
  zz_ov[i]=(zhist_big_ov[1][i+1]+zhist_big_ov[1][i])/2.
#we assume that the width of the intervals is always the same
dT=ate[1]-ate[0]
dm=am[1]-am[0]
dbeta=ab[1]-ab[0]
dntot=antot[1]-antot[0]
deovern=aeovern[1]-aeovern[0]
dsT3=asT3[1]-asT3[0]

nplots=len(pcount)
sp="    " #standard spacing in text output files
sfmt='{:8.4f}'
lfmt='{:12.7f}'

for p in range(nplots):
  if(pcount[p]>1):
    if(p<(nplots-1)):
      plabel=pnames[p]
      tlabel=psymb[p]
      fig, ((f1, f2), (f3, f4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6), dpi=100)   
      plt.subplots_adjust(hspace=0.6, wspace=0.3)
      fig.suptitle(tlabel)
      plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
      f1.set(title="time distribution",xlabel="t [fm/c]")
      f1.plot(tt,thist_cp[p][0]/np.sum(thist_cp[p][0]))
      f2.set(title="x distribution",xlabel="x [fm]")
      f2.plot(xx,xhist_cp[p][0]/np.sum(xhist_cp[p][0]))
      f3.set(title="y distribution",xlabel="y [fm]")
      f3.plot(yy,yhist_cp[p][0]/np.sum(yhist_cp[p][0]))
      f4.set(title="z distribution",xlabel="z [fm]")
      f4.plot(zz,zhist_cp[p][0]/np.sum(zhist_cp[p][0]))
      plt.tight_layout()
      plt.tight_layout()
      fig.savefig(od+"/check_plots_"+plabel+".png")
      plt.close()
      plabel=pnames[p]
      tlabel=psymb[p]
      fig, ((f1, f2), (f3, f4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6), dpi=100)   
      plt.subplots_adjust(hspace=0.6, wspace=0.3)
      fig.suptitle(tlabel+", |y|<"+str(raplim))
      plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
      f1.set(title="time distribution",xlabel="t [fm/c]")
      f1.plot(tt,thist[p][0]/np.sum(thist[p][0]))
      f2.set(title="x distribution",xlabel="x [fm]")
      f2.plot(xx,xhist[p][0]/np.sum(xhist[p][0]))
      f3.set(title="y distribution",xlabel="y [fm]")
      f3.plot(yy,yhist[p][0]/np.sum(yhist[p][0]))
      f4.set(title="z distribution",xlabel="z [fm]")
      f4.plot(zz,zhist[p][0]/np.sum(zhist[p][0]))
      plt.tight_layout()
      fig.savefig(od+"/check_plots_"+plabel+"_raplim.png")
      plt.close()
    else:
      plabel=pnames[p]
      tlabel=psymb[p]
      fig, ((f1, f2), (f3, f4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6), dpi=100)   
      plt.subplots_adjust(hspace=0.6, wspace=0.3)
      fig.suptitle(tlabel+", |y|<"+str(raplim))
      plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
      f1.set(title="time distribution",xlabel="t [fm/c]")
      f1.plot(tt,thist[p][0]/np.sum(thist[p][0]))
      f2.set(title="x distribution",xlabel="x [fm]")
      f2.plot(xx,xhist[p][0]/np.sum(xhist[p][0]))
      f3.set(title="y distribution",xlabel="y [fm]")
      f3.plot(yy,yhist[p][0]/np.sum(yhist[p][0]))
      f4.set(title="z distribution",xlabel="z [fm]")
      f4.plot(zz,zhist[p][0]/np.sum(zhist[p][0]))
      plt.tight_layout()
      fig.savefig(od+"/check_plots_"+plabel+"_raplim.png")
      plt.close()
      plabel=pnames[p]
      tlabel=psymb[p]
      fig, ((f1, f2), (f3, f4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6), dpi=100)   
      plt.subplots_adjust(hspace=0.6, wspace=0.3)
      plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
      fig.suptitle(tlabel)
      f1.set(title="time distribution",xlabel="t [fm/c]")
      f1.plot(tt,thist_cp[p][0]/np.sum(thist_cp[p][0]),ls="solid")
      f1.plot(tt,thist_big[0]/np.sum(thist_big[0]),ls="dashed")
      f2.set(title="x distribution",xlabel="x [fm]")
      f2.plot(xx,xhist_cp[p][0]/np.sum(xhist_cp[p][0]),ls="solid")
      f2.plot(xx,xhist_big[0]/np.sum(xhist_big[0]),ls="dashed")
      f3.set(title="y distribution",xlabel="y [fm]")
      f3.plot(yy,yhist_cp[p][0]/np.sum(yhist_cp[p][0]),ls="solid")
      f3.plot(yy,yhist_big[0]/np.sum(yhist_big[0]),ls="dashed")
      f4.set(title="z distribution",xlabel="z [fm]")
      f4.plot(zz,zhist_cp[p][0]/np.sum(zhist_cp[p][0]),ls="solid")
      f4.plot(zz,zhist_big[0]/np.sum(zhist_big[0]),ls="dashed")
      plt.tight_layout()
      fig.savefig(od+"/check_plots_"+plabel+".png")
      plt.close()
      plabel=pnames[p]
      tlabel=psymb[p]
      fig_ov, ((f1_ov, f2_ov), (f3_ov, f4_ov)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6), dpi=100)   
      plt.subplots_adjust(hspace=0.6, wspace=0.3)
      fig_ov.suptitle(tlabel)
      plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
      f1_ov.set(title="time distribution",xlabel="t [fm/c]")
      f1_ov.plot(tt,thist_cp[p][0]/np.sum(thist_cp[p][0]),ls="solid")
      f1_ov.plot(tt_ov,thist_big_ov[0]/np.sum(thist_big_ov[0]),ls="dashed")
      f2_ov.set(title="x distribution",xlabel="x [fm]")
      f2_ov.plot(xx,xhist_cp[p][0]/np.sum(xhist_cp[p][0]),ls="solid")
      f2_ov.plot(xx_ov,xhist_big_ov[0]/np.sum(xhist_big_ov[0]),ls="dashed")
      f3_ov.set(title="y distribution",xlabel="y [fm]")
      f3_ov.plot(yy,yhist_cp[p][0]/np.sum(yhist_cp[p][0]),ls="solid")
      f3_ov.plot(yy_ov,yhist_big_ov[0]/np.sum(yhist_big_ov[0]),ls="dashed")
      f4_ov.set(title="z distribution",xlabel="z [fm]")
      f4_ov.plot(zz,zhist_cp[p][0]/np.sum(zhist_cp[p][0]),ls="solid")
      f4_ov.plot(zz_ov,zhist_big_ov[0]/np.sum(zhist_big_ov[0]),ls="dashed")
      plt.tight_layout()
      fig_ov.savefig(od+"/check_plots_"+plabel+"_all_range.png")
      plt.close()
  

    NN=np.amax(phase[p][0])
    fig, ax= plt.subplots(dpi=400)
    plt.imshow(phase[p][0].transpose()/NN, extent=[phase[p][1][0], phase[p][1][-1], phase[p][2][0],phase[p][2][-1]], origin="lower",cmap='gnuplot2',aspect=3)
#    ax.set(title=tlabel)
    ax.set(xlabel='Baryon chemical potential [MeV]')
    ax.set(ylabel='Temperature [MeV]')
    plt.colorbar(label="Frozen hadrons (arbitrary units)",shrink=0.8)
    plt.savefig(od+"/T_mu_"+plabel+".png")
    plt.close()

    fig, ax = plt.subplots()
    ax.errorbar(arap,Trap[p][0],yerr=Trap[p][1],fmt='o')
#    ax.set(xlabel='y', ylabel='<T> [MeV]',  title=tlabel)
    ax.set(xlabel='y', ylabel='<T> [MeV]')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/T_rap_"+plabel+".png")
    plt.close()
    fp=open(od+"/T_rap_"+plabel+".dat","w")
    fp.write("#rapidity     temperature     std_dev_temperature\n")
    for wk in range(len(arap)):
        fp.write(sfmt.format(arap[wk])+sp+lfmt.format(Trap[p][0][wk])+sp+lfmt.format(Trap[p][1][wk])+"\n")
    fp.close()
    
    fig, ax = plt.subplots()
    ax.errorbar(arap,murap[p][0],yerr=murap[p][1],fmt='o')
#    ax.set(xlabel='y', ylabel='<$\mu_B$> [MeV]',  title=tlabel)
    ax.set(xlabel='y', ylabel='<$\mu_B$> [MeV]')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/mu_rap_"+plabel+".png")
    plt.close()
    fp=open(od+"/mu_rap_"+plabel+".dat","w")
    fp.write("#rapidity     mu_B     std_dev_mu_B\n")
    for wk in range(len(arap)):
        fp.write(sfmt.format(arap[wk])+sp+lfmt.format(murap[p][0][wk])+sp+lfmt.format(murap[p][1][wk])+"\n")
    fp.close()
    
    fig, ax = plt.subplots()
    ax.errorbar(arap,muoverTrap[p][0],yerr=muoverTrap[p][1],fmt='o')
#    ax.set(xlabel='y', ylabel='<$\mu_B$/T>',  title=tlabel)
    ax.set(xlabel='y', ylabel='<$\mu_B$/T>')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/muoverT_rap_"+plabel+".png")
    plt.close()
    fp=open(od+"/muoverT_rap_"+plabel+".dat","w")
    fp.write("#rapidity     mu/T     std_dev_mu/T\n")
    for wk in range(len(arap)):
        fp.write(sfmt.format(arap[wk])+sp+lfmt.format(muoverTrap[p][0][wk])+sp+lfmt.format(muoverTrap[p][1][wk])+"\n")
    fp.close()
    
    fig, ax = plt.subplots()
    ax.errorbar(apt,Tpt[p][0],yerr=Tpt[p][1],fmt='o')
#    ax.set(xlabel='$p_T$ [GeV]', ylabel='<T> [MeV]',  title=tlabel)
    ax.set(xlabel='$p_T$ [GeV]', ylabel='<T> [MeV]')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/T_pt_"+plabel+".png")
    plt.close()
    fp=open(od+"/T_pt_"+plabel+".dat","w")
    fp.write("#pT     temperature     std_dev_temperature\n")
    for wk in range(len(apt)):
        fp.write(sfmt.format(apt[wk])+sp+lfmt.format(Tpt[p][0][wk])+sp+lfmt.format(Tpt[p][1][wk])+"\n")
    fp.close()
    
    fig, ax = plt.subplots()
    ax.errorbar(apt,mupt[p][0],yerr=mupt[p][1],fmt='o')
#    ax.set(xlabel='$p_T$ [GeV]', ylabel='<$\mu_B$> [MeV]',  title=tlabel)
    ax.set(xlabel='$p_T$ [GeV]', ylabel='<$\mu_B$> [MeV]')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/mu_pt_"+plabel+".png")
    plt.close()
    fp=open(od+"/mu_pt_"+plabel+".dat","w")
    fp.write("#pT     mu_B     std_dev_mu_B\n")
    for wk in range(len(apt)):
        fp.write(sfmt.format(apt[wk])+sp+lfmt.format(mupt[p][0][wk])+sp+lfmt.format(mupt[p][1][wk])+"\n")
    fp.close()

    fig, ax = plt.subplots()
    ax.errorbar(apt,muoverTpt[p][0],yerr=muoverTpt[p][1],fmt='o')
#    ax.set(xlabel='$p_T$ [GeV]', ylabel='<$\mu_B$/T>',  title=tlabel)
    ax.set(xlabel='$p_T$ [GeV]', ylabel='<$\mu_B$/T>')
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/muoverT_pt_"+plabel+".png")
    plt.close()
    fp=open(od+"/muoverT_pt_"+plabel+".dat","w")
    fp.write("#pT     muoverTpt     std_dev_muoverTpt\n")
    for wk in range(len(apt)):
        fp.write(sfmt.format(apt[wk])+sp+lfmt.format(muoverTpt[p][0][wk])+sp+lfmt.format(muoverTpt[p][1][wk])+"\n")
    fp.close()

    NN=np.sum(temphist[p][0])
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel="T [MeV]",ylabel="1/N dN/dT [$MeV^{-1}$]")
    ax.set(xlabel="T [MeV]",ylabel="1/N dN/dT [$MeV^{-1}$]")
    ax.plot(ate,temphist[p][0]/(dT*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNdT_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNdT_"+plabel+".dat","w")
    fp.write("#T     1/N dN/dT     (N is "+str(NN)+")\n")
    for wk in range(len(ate)):
        fp.write(sfmt.format(ate[wk])+sp+lfmt.format(temphist[p][0][wk]/(dT*NN))+"\n")
    fp.close()

    NN=np.sum(muhist[p][0])
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel="$\mu_B$ [MeV]",ylabel="1/N dN/d$\mu_B$ [$MeV^{-1}$]")
    ax.set(xlabel="$\mu_B$ [MeV]",ylabel="1/N dN/d$\mu_B$ [$MeV^{-1}$]")
    ax.plot(am,muhist[p][0]/(dm*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNdmu_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNdmu_"+plabel+".dat","w")
    fp.write("#mu_B     1/N dN/dmu_B     (N is "+str(NN)+")\n")
    for wk in range(len(am)):
        fp.write(sfmt.format(am[wk])+sp+lfmt.format(muhist[p][0][wk]/(dm*NN))+"\n")
    fp.close()

    NN=np.sum(ntot[p][0])
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel="$\mu_B$ [MeV]",ylabel="1/N dN/d$\mu_B$ [$MeV^{-1}$]")
    ax.set(xlabel=r"$\rho_{b+ab}$ [1/fm^3]",ylabel="1/N dN/d"+r"$\rho_{b+ab}$ [$fm^3$]")
    ax.plot(antot,ntot[p][0]/(dntot*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNdntot_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNdntot_"+plabel+".dat","w")
    fp.write("#n_b+n_ab     1/N dN/d(n_b+n_ab)     (N is "+str(NN)+")\n")
    for wk in range(len(antot)):
        fp.write(sfmt.format(antot[wk])+sp+lfmt.format(ntot[p][0][wk]/(dntot*NN))+"\n")
    fp.close()

    NN=np.sum(e_over_n[p][0])
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel="$\mu_B$ [MeV]",ylabel="1/N dN/d$\mu_B$ [$MeV^{-1}$]")
    ax.set(xlabel="e/n [GeV]",ylabel="1/N dN/d(e/n) [$GeV^{-1}$]")
    ax.plot(aeovern,e_over_n[p][0]/(deovern*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNde_over_n_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNde_over_n_"+plabel+".dat","w")
    fp.write("#n/e     1/N dN/d(n/e)     (N is "+str(NN)+")\n")
    for wk in range(len(aeovern)):
        fp.write(sfmt.format(aeovern[wk])+sp+lfmt.format(e_over_n[p][0][wk]/(deovern*NN))+"\n")
    fp.close()

    NN=np.sum(sT3[p][0])
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel="$\mu_B$ [MeV]",ylabel="1/N dN/d$\mu_B$ [$MeV^{-1}$]")
    ax.set(xlabel="$s/T^3$",ylabel="1/N dN/d$(s/T^3)$")
    ax.plot(asT3,sT3[p][0]/(dsT3*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNdsT3_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNdsT3_"+plabel+".dat","w")
    fp.write("#s/T^3     1/N dN/d(s/T^3)    (N is "+str(NN)+")\n")
    for wk in range(len(antot)):
        fp.write(sfmt.format(asT3[wk])+sp+lfmt.format(sT3[p][0][wk]/(dsT3*NN))+"\n")
    fp.close()
  
    NN=np.sum(betahist[p][0]) 
    fig, ax = plt.subplots()
#    ax.set(title=tlabel,xlabel=r"$\beta$ [c units]",ylabel="1/N dN/d"+r"$\beta$")
    ax.set(xlabel=r"$\beta$ [c units]",ylabel="1/N dN/d"+r"$\beta$")
    ax.plot(ab,betahist[p][0]/(dbeta*NN))
    plt.grid(b=True, which='major', color='gainsboro', linestyle='-')
    fig.savefig(od+"/dNdbeta_"+plabel+".png")
    plt.close()
    fp=open(od+"/dNdbeta_"+plabel+".dat","w")
    fp.write("#beta_transverse     1/N dN/dbeta     (N is "+str(NN)+")\n")
    for wk in range(len(ab)):
        fp.write(sfmt.format(ab[wk])+sp+lfmt.format(betahist[p][0][wk]/(dbeta*NN))+"\n")
    fp.close()
