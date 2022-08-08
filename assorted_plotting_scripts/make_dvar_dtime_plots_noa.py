# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3
import numpy as np
import matplotlib.pyplot as plt
import pickle


files=("timedis_elb_1_23noa_data.pickle","timedis_elb_4noa_data.pickle","timedis_elb_10_8noa_data.pickle","timedis_ecm_7_7_data.pickle","timedis_ecm_9_1_data.pickle","timedis_ecm_11_5_data.pickle","timedis_ecm_14_5_data.pickle","timedis_ecm_19_6_data.pickle", "timedis_ecm_27_data.pickle","timedis_ecm_39_data.pickle","timedis_ecm_62_4_data.pickle")
colors=("xkcd:neon yellow","lime", "springgreen", "c", "b", "xkcd:blue purple","m", "xkcd:shocking pink", "r", "brown", "k")
labels=(r'$E_{lab}$=1.23 AGeV',r'$E_{lab}$=4.0 AGeV',r'$E_{lab}$=10.8 AGeV',r'$\sqrt{s_{NN}}$=7.7 GeV',r'$\sqrt{s_{NN}}$=9.1 GeV',r'$\sqrt{s_{NN}}$=11.5 GeV',r'$\sqrt{s_{NN}}$=14.5 GeV',r'$\sqrt{s_{NN}}$=19.6 GeV',r'$\sqrt{s_{NN}}$=27 GeV',r'$\sqrt{s_{NN}}$=39 GeV',r'$\sqrt{s_{NN}}$=62.4 GeV')

s_list=[]
dsdt_list=[]
nbab_list=[]
eovn_list=[]
sT3_list=[]

#we load the files individually
for k in range(len(files)):
    with open(files[k],"rb") as pi:
         #we assume that all the histograms use the same bins, so we don't care to save separately the grid informations
         xsize,ysize,zsize,dt,tmax,nt,s_res,dsdt_res,eovn_res,nbab_res,sT3_res=pickle.load(pi)[:]
         s_list.append(s_res.copy())
         dsdt_list.append(dsdt_res.copy())
         nbab_list.append(nbab_res.copy())
         eovn_list.append(eovn_res.copy())
         sT3_list.append(sT3_res.copy())

data=None #we free the memory occupied by data
tarr=np.linspace(dt/2.,tmax-dt/2.,num=nt,endpoint=True)

#plot of s(t)
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,8.2,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|x|,|y|,|z|<5 fm') 
ax.set_xlabel('t [fm]', fontsize=14)
ax.set_ylabel('s [$fm^{-3}$]', fontsize=14)
ax.set_xlim(0.,tmax)
ax.set_ylim(0.,20)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(tarr,s_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("s_vs_time_noa.png",dpi=300)
plt.close()

#plot of ds/dt(t)
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,2.2,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|x|,|y|,|z|<5 fm') 
ax.set_xlabel('t [fm]', fontsize=14)
ax.set_ylabel('ds/dt [$fm^{-4}$]', fontsize=14)
ax.set_xlim(0.,tmax)
ax.set_ylim(-6,4)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(tarr,dsdt_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("dsdt_vs_time_noa.png",dpi=300)
plt.close()

#plot of n_bab(t)
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,0.8,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|x|,|y|,|z|<5 fm') 
ax.set_xlabel('t [fm]', fontsize=14)
ax.set_ylabel(r'$\rho_{b+ab}$'+' [$fm^{-3}$]', fontsize=14)
ax.set_xlim(0.,tmax)
ax.set_ylim(0.,1.6)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(tarr,nbab_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("nbab_vs_time_noa.png",dpi=300)
plt.close()

#plot of e/n(t)
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,0.2,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|x|,|y|,|z|<5 fm') 
ax.set_xlabel('t [fm]', fontsize=14)
ax.set_ylabel('<E>/<N> [GeV]', fontsize=14)
ax.set_xlim(0.,tmax)
ax.set_ylim(0.,2.25)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(tarr,eovn_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("energy_per_particle_vs_time_noa.png",dpi=300)
plt.close()

#plot of s/T^3
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,16,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|x|,|y|,|z|<5 fm') 
ax.set_xlabel('t [fm]', fontsize=14)
ax.set_ylabel('$s/T^3$', fontsize=14)
ax.set_xlim(0.,tmax)
ax.set_ylim(0.,25)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(tarr,sT3_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("sT3_vs_time_noa.png",dpi=300)
plt.close()
