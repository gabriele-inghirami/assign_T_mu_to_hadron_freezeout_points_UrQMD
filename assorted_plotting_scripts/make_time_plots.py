# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3
# version 0.5 - 27/05/2020

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip


#distance from the center of the grid
xsize=5
ysize=5
zsize=5

#time resolution
dt=0.2

#time max
tmax=30

#number of timesteps (automatically set from dt and tmax)
nt=int(math.floor(tmax/dt))

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<2):
  print ('Syntax: ./get_info.py <inputfile> <outputlabel>')
  sys.exit(1)

inputfile=sys.argv[1]
outputlabel=sys.argv[2]

#first, we check that files exist
if(not(os.path.isfile(inputfile))):
  print(inputfile+" does not exist. I quit.\n")
  sys.exit(1)

s_arr=np.zeros((nt,2),dtype=np.float64)
e_over_n_arr=np.zeros((nt,2),dtype=np.float64)
nbab_arr=np.zeros((nt,2),dtype=np.float64)
sT3_arr=np.zeros((nt,2),dtype=np.float64)

s_res=np.zeros(nt,dtype=np.float64)
dsdt_res=np.zeros(nt,dtype=np.float64)
eovn_res=np.zeros(nt,dtype=np.float64)
nbab_res=np.zeros(nt,dtype=np.float64)
sT3_res=np.zeros(nt,dtype=np.float64)
dsdt_res=np.zeros(nt,dtype=np.float64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")

hc3=197.326**3
lines=0
hadrons=0
for line in pi:
    #decomment the following line to test only a limited subset of the input data
    #if(lines==10000):
    #   break
    stuff=line.split() 
    lines=lines+1
    t,x,y,z=np.float64(stuff[3:7])
    if(( t < tmax) and (abs(x) < xsize) and (abs(y) < ysize) and (abs(z)< zsize)):
      hadrons=hadrons+1
      i=int(math.floor(t/dt))
      temp=np.float(stuff[10])
      en,num_had,s,rho_bab=np.float64(stuff[12:16])
      s_arr[i,0]=s_arr[i,0]+s
      s_arr[i,1]=s_arr[i,1]+1

      e_over_n_arr[i,0]=e_over_n_arr[i,0]+en
      e_over_n_arr[i,1]=e_over_n_arr[i,1]+num_had
  
      nbab_arr[i,0]=nbab_arr[i,0]+rho_bab
      nbab_arr[i,1]=nbab_arr[i,1]+1

      if(temp !=0 ):
        sT3_arr[i,0]=sT3_arr[i,0]+s/temp**3 #later we multiply also the whole array by hbarc^3
        sT3_arr[i,1]=sT3_arr[i,1]+1

pi.close()

sT3_arr[:,0]=sT3_arr[:,0]*hc3

#now we print the results into a file
sp="          "
fout=open(outputlabel+"_data.txt","w")
fout.write("#Limits: t < "+'{:5.2f}'.format(t)+"fm, |x|<"+'{:5.2f}'.format(x)+"fm, |y|<"+'{:5.2f}'.format(y)+"fm, |y|<"+'{:5.2f}'.format(y)+"\n")
fout.write("#t       <s>(fm^-3)       ds/dt(fm^-4)       <E>/<N>(GeV)     <n>_{b+ab}(fm^-3)     <s/T^3>\n")
for i in range(0,nt):
    t=(i+0.5)*dt
    if(s_arr[i,1]==0):
      s=0.
    else:
      s=s_arr[i,0]/s_arr[i,1]
    s_res[i]=s
    if((i==0) or (s_arr[i,1]==0) or (s_arr[i-1,1]==0)): #the dsdt array has one cell less than the others
      dsdt=0
    else:
      dsdt=(s_arr[i,0]/s_arr[i,1]-s_arr[i-1,0]/s_arr[i-1,1])/dt
    dsdt_res[i]=dsdt
    if(e_over_n_arr[i,1]==0):
      e_over_n=0.
    else:
      e_over_n=e_over_n_arr[i,0]/e_over_n_arr[i,1]
    eovn_res[i]=e_over_n
    if(nbab_arr[i,1]==0):
      nbab=0.
    else:
      nbab=nbab_arr[i,0]/nbab_arr[i,1]
    nbab_res[i]=nbab
    if(sT3_arr[i,1]==0):
      sT3=0.
    else:
      sT3=sT3_arr[i,0]/sT3_arr[i,1]
    sT3_res[i]=sT3

    fout.write('{:5.2f}'.format(t)+sp+'{:7.4f}'.format(s)+sp+'{:7.4f}'.format(dsdt)+sp+'{:7.4f}'.format(e_over_n)+sp+'{:7.4f}'.format(nbab)+sp+'{:7.4f}'.format(sT3)+"\n")
    
fout.close()

with open(outputlabel+"_data.pickle","wb") as po:
     pickle.dump((xsize,ysize,zsize,dt,tmax,nt,s_res,dsdt_res,eovn_res,nbab_res,sT3_res),po)
