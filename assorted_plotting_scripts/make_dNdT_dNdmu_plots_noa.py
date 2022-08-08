# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3

import numpy as np
import matplotlib.pyplot as plt
import pickle


files=("res_elb_1_23_Fnoa_plot_data_binary.pickle","res_elb_4_Fnoa_plot_data_binary.pickle","res_elb_10_8_Fnoa_plot_data_binary.pickle","res_ecm_7_7_Fnoa_plot_data_binary.pickle","res_ecm_9_1_Fnoa_plot_data_binary.pickle","res_ecm_11_5_Fnoa_plot_data_binary.pickle","res_ecm_14_5_Fnoa_plot_data_binary.pickle","res_ecm_19_6_Fnoa_plot_data_binary.pickle", "res_ecm_27_Fnoa_plot_data_binary.pickle","res_ecm_39_Fnoa_plot_data_binary.pickle","res_ecm_62_4_Fnoa_plot_data_binary.pickle")
colors=("xkcd:neon yellow","lime", "springgreen", "c", "b", "xkcd:blue purple","m", "xkcd:shocking pink", "r", "brown", "k")
labels=(r'$E_{lab}$=1.23 AGeV',r'$E_{lab}$=4.0 AGeV',r'$E_{lab}$=10.8 AGeV',r'$\sqrt{s_{NN}}$=7.7 GeV',r'$\sqrt{s_{NN}}$=9.1 GeV',r'$\sqrt{s_{NN}}$=11.5 GeV',r'$\sqrt{s_{NN}}$=14.5 GeV',r'$\sqrt{s_{NN}}$=19.6 GeV',r'$\sqrt{s_{NN}}$=27 GeV',r'$\sqrt{s_{NN}}$=39 GeV',r'$\sqrt{s_{NN}}$=62.4 GeV')


dNdT_list=[]
dNdmu_list=[]

#we load the files individually
#the index 25 and 26 corresponds to the histograms of dN/dT and dN/dmu, respectively. We take the last index [-1], corresponding to all particles, and then 0 and 1 to get the counts and the border of the bins
for k in range(len(files)):
    with open(files[k],"rb") as pi:
         data=pickle.load(pi)
         if(k==0):
           #we compute the centers of the bins and we save them into the arrays temp and mu
           #we assume that all the histograms use the same bins, so we do this step only once
           lte=len(data[25][-1][1])-1
           temp=np.zeros(lte)
           lm=len(data[26][-1][1])-1
           mu=np.zeros(lm)
           for i in range(lte):
               temp[i]=(data[25][-1][1][i+1]+data[25][-1][1][i])/2.
           for i in range(lm):
               mu[i]=(data[26][-1][1][i+1]+data[26][-1][1][i])/2.
           dT=temp[1]-temp[0]
           dmu=mu[1]-mu[0]
         dNdT_list.append(data[25][-1][0].copy()/(dT*np.sum(data[25][-1][0])))
         dNdmu_list.append(data[26][-1][0].copy()/(dmu*np.sum(data[26][-1][0])))


data=None #we free the memory occupied by data

#plot of 1/N dN/dT
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(10,0.035,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|y|$\leq$0.5\np$_\mathrm{T}\leq$3 GeV') 
ax.set_xlabel('T [MeV]', fontsize=14)
ax.set_ylabel('1/N dN/dT [MeV$^{-1}$]', fontsize=14)
ax.set_xlim(0.,260.)
ax.set_ylim(0.,0.05)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(temp,dNdT_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("dN_over_NdT_noa.png",dpi=300)
plt.close()

#plot of 1/N dN/dmu
fig, ax=plt.subplots(figsize=(6.5,6.5))
ax.set_title(r'$\pi$', fontsize=14)
ax.grid(b=True, which='major', color='grey', linestyle='-', alpha=0.5)
ax.grid(b=True, which='minor', color='grey', linestyle=':', alpha=0.5)
ax.text(400,0.0205,'Au+Au (UrQMD)\nb$\leq$3.4 fm\n|y|$\leq$0.5\np$_\mathrm{T}\leq$3 GeV') 
ax.set_xlabel(r'$\mu$ [MeV]', fontsize=14)
ax.set_ylabel(r'1/N dN/d$\mu$ [MeV$^{-1}$]', fontsize=14)
ax.set_xlim(0.,1100.)
ax.set_ylim(0.,0.030)
#ax.set_yscale('log')
for i in range(len(files)):
    ax.plot(mu,dNdmu_list[i], color=colors[i],label=labels[i])
ax.legend()
fig.savefig("dN_over_Ndmu_noa.png",dpi=300)
plt.close()
