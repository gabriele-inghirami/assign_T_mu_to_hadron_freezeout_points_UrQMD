# version 2.2.0 - 25/05/2020

# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3

#this version works with the files generated by store_cgnew_v2.1.1.py, in which the pressure, hadron, entropy and total baryon+antibaryon densities are saved, too
#in this version the energy density that it is printed in the particle list is the same used in the EoS (standard or the anisotropy corrected)
#this version creates also histograms with respect to nbab (total baryon + antibaryon density), e/n (energy per particle) and s/T^3

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip


"""It analyzes the pickle binary data produced by store_cg.py (v. 2.1), also if gzipped."""

#set true to compute only T and mu without computing also the plots
onlydata=False

#it chooses whether to use standard (False) or the anisotropy corrected (True) energy density
use_aniso=False

#it tells how many columns have the lines with data (4 or 5)
datacolumns=5 

#it sets a limit on the ratio between transverse and longitudinal pressure: abs(log10(P_parallel/P_transverse))<tranverse_parallel_pressure_limit
#use a 0 or negative value to disable the limit
tranverse_parallel_pressure_limit=-0.6

#box size for plots (it can be larger than the coarse graining box, but it would be useless)
#for each of the spatial dimensions, the box is 2 times larger, here we consider only x,y,z>0
xsize=40
ysize=40
zsize=30
tsize=40

#bins for the histograms
xcells=80
ycells=80
zcells=60
tcells=160
#number of bins and ranges multiplying factor for the overall histograms
ovf=1

#transverse momentum (units: GeV) and rapidity ranges and bins
ptmax=3
ptlim=ptmax
ptbins=20
rapmax=4 #rapmax is the maximum accepted rapidity, used also for rapidity plots
rapbins=20
raplim=0.2 #in the plots at central rapidity, |y|<raplim

#transverse velocity (units: c) bins
beta_bins=40 

#for the plots (units: MeV)
min_baryon_chem_to_plot=0
max_baryon_chem_to_plot=1000
baryon_chem_bins=125
min_temperature_to_plot=0
max_temperature_to_plot=350
temperature_bins=70

max_nbab_to_plot=0.5
nbab_bins=100
max_e_over_n=1.5
e_over_n_bins=75
max_sT3=25
sT3_bins=100

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<3):
   print ('Syntax: ./get_info.py <inputfile> <outputfile> <coarse_file_data 1> [coarse_file_data 2] ...')
   sys.exit(1)

inputfile=sys.argv[1]
outputfile=sys.argv[2]
coarsefiles=sys.argv[3:]

maxparticles=0 #maximum number of particle to place on the coarse grained grid
with open(inputfile) as f:
    i=-1
    for i, l in enumerate(f):
        pass
maxparticles=i + 1

#NOTE: if you add another index, please, remind to increase the number of elements in the partdata_big, partdata and partdata_det below
# Moreover, please, pay attention to the last index copied from partdata_big to partdata (at the moment, kbbeta)
#we define indexes for convenience in locating the data in partdata_big, partdata and partdata_raplim
kbid=0;kbt=1;kbx=2;kby=3;kbz=4;kbvx=5;kbvy=6;kbvz=7;kbmu=8;kbtemp=9;kbrho=10;kben=11;kbrap=12;kbpt=13;kbbeta=14;kbnbab=15;kbeovern=16;kbsT3=17;kbpnum=18;kbpch=19;kbparent=20;
elements_of_partdata_big=21
#now we define indexes for convenience in locating the data in partdata_det and partdata_det_raplim
kt=0;kx=1;ky=2;kz=3;kvx=4;kvy=5;kvz=6;kmu=7;ktemp=8;krho=9;ken=10;krap=11;kpt=12;kbeta=13;knbab=14;keovern=15;ksT3=16;
elements_of_partdata_det=17
elements_of_partdata=elements_of_partdata_det+1 

ptype={(1,0):(0,"N",0.939565413),(-1,0):(1,"aN",0.939565413),(1,1):(2,"p",0.938272081),(-1,-1):(3,"ap",0.938272081),(27,0):(4,"L",1.115683),(-27,0):(5,"aL",1.115683),(40,0):(6,"S0",1.192642),(-40,0):(7,"aS0",1.192642),(49,0):(8,"Xi0",1.31486),(-49,0):(9,"aXi0",1.31486),(55,-1):(10,"Om",1.67245),(-55,+1):(11,"aOm",1.67245),(101,1):(12,"pip",0.13957061),(101,0):(13,"pi0",0.1349770),(101,-1):(14,"pim",0.13957061),(106,1):(15,"Kp",0.493677),(-106,-1):(16,"Km",0.493677),(106,0):(17,"K0",0.497611),(-106,0):(18,"aK0",0.497611),(99,1):(19,"deu",1.875613),(-99,-1):(20,"adeu",1.875613)}

ntype=len(ptype)

#arrays to store coarse graining data
print("Allocating and initializing arrays... ")
partdata_big=np.zeros((maxparticles,elements_of_partdata_big))
pcount=np.zeros(ntype+1,dtype=np.int)
pcount_raplim=np.zeros(ntype+1,dtype=np.int)
print("Done.")

sp="      "
fp='{:>8.4f}'
fn='{:>9.6f}'
iinum='{:>5d}'
hc3=197.326**3

#we define an auxiliary function to idetify integers
def test_integer(xyz):
     try:
        bitbucket=int(xyz)
     except ValueError:
        return False
     return True
 
particles_out_of_box=0
particles_out_of_time=0
particles_out_of_box_raplim=0
particles_out_of_time_raplim=0
tot_particles=0
tot_particles_raplim=0
empty_cells=0
empty_cells_raplim=0
aniso_exceeding_cells=0
good=0
good_raplim=0
unknown=0
inbox=0
inbox_raplim=0
intime=0
intime_raplim=0

#now we read the position of the deutoreon and we assign the corresponding temperature and mu
indata = open(inputfile,"r")

outdata = open(outputfile+".dat","w")

#out_unknown=open(outputfile+"_unknown_part_list.dat","w")
#out_unknown.write("Particle of unknown types:\nLine  id  charge\n")

count_lines=0
for line in indata:
    count_lines=count_lines+1
    stuff_raw=line.split()
    ncolumns=len(stuff_raw)
    if(ncolumns==datacolumns):
      try:
         if(datacolumns==4):
           t,x,y,z=np.float64(stuff_raw[:])
         else:
           t,x,y,z,bitbucket=np.float64(stuff_raw[:])
      except(ValueError):
         print("Reading error at line: "+str(count_lines))
         print(line)
         continue
      if(not(np.isfinite(t) and np.isfinite(x) and np.isfinite(y) and np.isfinite(z))):
         continue
      p0=0.
      px=0.
      py=0.
      pz=0.
      part_rap=0.
      part_pt=0.      
      pid=14
      if((abs(part_rap)<raplim) and (part_pt<ptmax)):
        tot_particles_raplim=tot_particles_raplim+1

      partdata_big[tot_particles,kbid]=pid
      partdata_big[tot_particles,kbt]=t
      partdata_big[tot_particles,kbx]=x
      partdata_big[tot_particles,kby]=y
      partdata_big[tot_particles,kbz]=z
      partdata_big[tot_particles,kbpt]=part_pt
      partdata_big[tot_particles,kbrap]=part_rap
      partdata_big[tot_particles,kbpnum]=np.float64(stuff_raw[0])
      partdata_big[tot_particles,kbpch]=np.float64(stuff_raw[1])
      partdata_big[tot_particles,kbparent]=np.float64(stuff_raw[2])
      tot_particles=tot_particles+1
      if(tot_particles%500000 == 0):
         print("Hadrons read so far: "+str(tot_particles))

indata.close()
#out_unknown.close()

print("Resizing particle array.")
partdata_tmp=partdata_big[0:tot_particles,:]
partdata_big=None
partdata_big=partdata_tmp
partdata_tmp=None

ncofiles=len(coarsefiles)
if(ncofiles==1):
  print("Reading coarse data file... ")
else:
  print("Reading coarse data files... ")

first_time=True

for df in range(ncofiles):
  if(coarsefiles[df][-3:]==".gz"):
    print("Opening gzipped file "+coarsefiles[df])
    pi=gzip.open(coarsefiles[df],"rb")
  else:
    print("Opening file "+coarsefiles[df])
    pi=open(coarsefiles[df],"rb")
     
  data=pickle.load(pi)

  intt,inxx,inyy,inzz,invx,invy,invz,inrho,inen,inmuSTD,intempSTD,inmuANI,intempANI,ptra,ppar,num_had,pressSTD,sSTD,pressANI,sANI,rho_bab=data[:]

  tmax=intt[-1]
  tmin=intt[0]  
  print("Tmin and tmax are: "+str(tmin)+"  "+str(tmax))
  dt=intt[1]-intt[0]
  dx=inxx[1]-inxx[0]
  dy=inyy[1]-inyy[0]
  dz=inzz[1]-inzz[0]

  xmin=inxx[0]-0.5*dx
  ymin=inyy[0]-0.5*dy
  zmin=inzz[0]-0.5*dz

  nx=len(inxx)
  ny=len(inyy)
  nz=len(inzz)
  
  if(use_aniso):
    print("Use temperature and baryon chemical potential with anisotropic correction")
    inmu=inmuANI
    intemp=intempANI
  else:
    print("Use standard temperature and baryon chemical potential, without anisotropic correction")
    inmu=inmuSTD
    intemp=intempSTD

  if(tranverse_parallel_pressure_limit>0):
    aniso_limit=True
  else:
    aniso_limit=False

  #cdata.close()
  print("File "+coarsefiles[df]+" read, now I am associating T and mu to the particles in the computational box included in the file...")

  for pp in range(0,tot_particles):
    pid=int(partdata_big[pp,kbid])
    t=partdata_big[pp,kbt]
    x=partdata_big[pp,kbx]
    y=partdata_big[pp,kby]
    z=partdata_big[pp,kbz]
    part_pt=partdata_big[pp,kbpt]
    part_rap=partdata_big[pp,kbrap]
    if((t>=tmin-dt/2.) and (t<tmax+dt/2.)):
      intime=intime+1
      if((abs(part_rap)<raplim) and (part_pt<ptlim)):
        intime_raplim=intime_raplim+1
      h=int(math.floor((t+dt/2.-tmin)/dt))
      i=int(math.floor((x-xmin)/dx))
      if((i >= 0) and (i<nx)):
        j=int(math.floor((y-ymin)/dy))
        if((j >= 0) and (j<ny)):
          k=int(math.floor((z-zmin)/dz))
          if((k >= 0) and (k<nz)):
            inbox=inbox+1
            if((abs(part_rap)<raplim) and (part_pt<ptlim)):
              inbox_raplim=inbox_raplim+1
            if(aniso_limit and ((ppar[h,i,j,k]==0) or (ptra[h,i,j,k]==0))):
              skip_because_aniso=True
            elif(aniso_limit and (abs(np.log10(ppar[h,i,j,k]/ptra[h,i,j,k]))>=tranverse_parallel_pressure_limit)):
              aniso_exceeding_cells=aniso_exceeding_cells+1
              skip_because_aniso=True
            else:
              skip_because_aniso=False
            if((intemp[h,i,j,k]==0) or (inmu[h,i,j,k]==0) or skip_because_aniso):
              empty_cells=empty_cells+1
              if(abs(part_rap)<raplim):
                 empty_cells_raplim=empty_cells_raplim+1
            else:
              if(use_aniso): #in this case we renormalize the energy density, as that one that we read is always the "standard" one
                aniso_ratio=ppar[h,i,j,k]/ptra[h,i,j,k] #if we are here skip_because_aniso is False and therefore ptra!=0
                x_aniso=aniso_ratio**(4./3.)
                if(x_aniso<1):
                   r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctanh(np.sqrt(1.-x_aniso)))/(np.sqrt(1.-x_aniso)))
                elif(x_aniso>1):
                   r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctan(np.sqrt(x_aniso-1.)))/(np.sqrt(x_aniso-1.)))
                else: #x_aniso==1
                   r_aniso=1.
                energy_density=inen[h,i,j,k]/r_aniso
                entropy_density=sANI[h,i,j,k]
              else:
                energy_density=inen[h,i,j,k]
                entropy_density=sSTD[h,i,j,k]

              partdata_big[pp,kbvx]=invx[h,i,j,k]
              partdata_big[pp,kbvy]=invy[h,i,j,k]
              partdata_big[pp,kbvz]=invz[h,i,j,k]
              partdata_big[pp,kbbeta]=math.sqrt(invx[h,i,j,k]**2+invy[h,i,j,k]**2)
              partdata_big[pp,kbmu]=inmu[h,i,j,k]*1000
              partdata_big[pp,kbtemp]=intemp[h,i,j,k]*1000
              partdata_big[pp,kbrho]=inrho[h,i,j,k]
              partdata_big[pp,kben]=energy_density
              partdata_big[pp,kbnbab]=rho_bab[h,i,j,k]
              if(num_had[h,i,j,k]!=0):
                 partdata_big[pp,kbeovern]=energy_density/num_had[h,i,j,k]
                 #print("en dens, num_had and e/n: "+str(energy_density)+"  "+str(num_had[h,i,j,k])+"  "+str(partdata_big[pp,kbeovern]))
              if(entropy_density!=0):
                 partdata_big[pp,kbsT3]=entropy_density*hc3/(partdata_big[pp,kbtemp])**3
                 #print("entr dens, T and sT3: "+str(entropy_density)+"  "+str(partdata_big[pp,kbtemp])+"  "+str(partdata_big[pp,kbsT3])+"  "+str(hc3))
                
              if((pid>=0) and (abs(part_rap)<rapmax) and (part_pt<ptlim)):
#                print(str(pid))
                pcount[pid]=pcount[pid]+1
                good=good+1
                if(abs(part_rap)<raplim):
                  pcount_raplim[pid]=pcount_raplim[pid]+1
                  good_raplim=good_raplim+1
              outdata.write(iinum.format(int(partdata_big[pp,kbpnum]))+sp+iinum.format(int(partdata_big[pp,kbpch]))+sp+iinum.format(int(partdata_big[pp,kbparent]))+sp+fn.format(t)+sp+fn.format(x)+sp+fn.format(y)+sp+fn.format(z)+sp+fn.format(part_rap)+sp+fn.format(part_pt)+sp+fn.format(inmu[h,i,j,k]*1000)+sp+fn.format(intemp[h,i,j,k]*1000)+sp+fn.format(inrho[h,i,j,k])+sp+fn.format(energy_density)+sp+fn.format(num_had[h,i,j,k])+sp+fn.format(entropy_density)+sp+fn.format(rho_bab[h,i,j,k])+"\n")
                          
  pi.close()

outdata.close()

particles_out_of_box=tot_particles-inbox
particles_out_of_time=tot_particles-intime
particles_out_of_box_raplim=tot_particles_raplim-inbox_raplim
particles_out_of_time_raplim=tot_particles_raplim-intime_raplim

print("Done.")
print("Total particles: "+str(tot_particles))
print("Particles of unknown type: "+str(unknown))
print("Particles out of the space-time box: "+str(particles_out_of_box)+", ("+str(float(particles_out_of_box/tot_particles)*100)+"%)")
print("Particles out of the time box: "+str(particles_out_of_time)+", ("+str(float(particles_out_of_time/tot_particles)*100)+"%)")
print("Particles out of the space-time box with |y|<"+str(raplim)+": "+str(particles_out_of_box_raplim)+", ("+str(float(particles_out_of_box_raplim/tot_particles_raplim)*100)+"% of all particles with |y|<"+str(raplim)+")")
print("Particles out of the time box with |y|<"+str(raplim)+" and pT<"+str(ptlim)+": "+str(particles_out_of_time_raplim)+", ("+str(float(particles_out_of_time_raplim/tot_particles_raplim)*100)+"% of all particles with |y|<"+str(raplim)+")")
print("Particles in the space-time box: "+str(inbox))
print("Particles in the space-time box with |y|<"+str(raplim)+" and pT<"+str(ptlim)+": "+str(inbox_raplim))
print("Particles in empty coarse-graining cells: "+str(empty_cells)+", ("+str(float(empty_cells/inbox)*100)+"%)")
print("Particles in empty coarse-graining cells with |y|<"+str(raplim)+" and pT<"+str(ptlim)+": "+str(empty_cells_raplim)+", ("+str((float(empty_cells_raplim/inbox_raplim))*100)+"% of all particles with |y|<"+str(raplim))
if(aniso_limit):
  print("Particle discarded because they exceed the anisotropy limit ("+str(tranverse_parallel_pressure_limit)+"): "+str(aniso_exceeding_cells))
print("Particles with associated T, mu: "+str(good)+" , with |y|<"+str(raplim)+" and pT<"+str(ptlim)+": "+str(good_raplim))
partdata_det=[]
partdata_det_raplim=[]
for k, v in ptype.items():
    numpart=pcount[v[0]] 
    numpart_raplim=pcount_raplim[v[0]]
    print(v[1]+" : "+str(numpart)+", with |y|<"+str(raplim)+": "+str(numpart_raplim))
    partdata_det.append(np.zeros((numpart,elements_of_partdata_det),dtype=np.float64))
    partdata_det_raplim.append(np.zeros((numpart_raplim,elements_of_partdata_det),dtype=np.float64))

print("\n")
if(onlydata):
  print("You chose to compute only T and mu for the particle data file, not the plots. I already did it and I stop here.")
  sys.exit(0)

print("Creating array with ''good'' particles")
partdata=np.zeros((good,elements_of_partdata),dtype=np.float64)
partdata_raplim=np.zeros((good_raplim,elements_of_partdata),dtype=np.float64)
j=0
k=0
for i in range(0,tot_particles):
    if((partdata_big[i,kbtemp]>0) and (partdata_big[i,kbmu]>0) and (partdata_big[i,kbid]>=0) and (abs(partdata_big[i,kbrap])<rapmax) and (partdata_big[i,kbpt]<ptlim)):
      partdata[j,:]=partdata_big[i,0:elements_of_partdata]
      j=j+1
      if(abs(partdata_big[i,kbrap])<raplim):
        partdata_raplim[k,:]=partdata_big[i,0:elements_of_partdata]
        k=k+1

print("Reallocating arrays for histograms.")
intt=None
inxxx=None
inyy=None
inzz=None
invx=None
invy=None
invz=None
inmu=None
intemp=None
inrho=None
inen=None
inmuSTD=None
intempSTD=None
inmuANI=None
intempANI=None
ptra=None
ppar=None
num_had=None
pressSTD=None
sSTD=None
pressANI=None
sANI=None
rho_bab=None
pcount[-1]=good

print("Creating bins for histograms")
#we create the bins. We will compute histograms with numpy, therefore we need to provide the edges of the bins
tt=np.linspace(tmin-dt/2.,tsize+dt/2,num=tcells+1,endpoint=True)
xx=np.linspace(-xsize,xsize,num=xcells+1,endpoint=True)
yy=np.linspace(-ysize,ysize,num=ycells+1,endpoint=True)
zz=np.linspace(-zsize,zsize,num=zcells+1,endpoint=True)
mb=np.linspace(min_baryon_chem_to_plot,max_baryon_chem_to_plot,num=baryon_chem_bins+1,endpoint=True)
tb=np.linspace(min_temperature_to_plot,max_temperature_to_plot,num=temperature_bins+1,endpoint=True)
bb=np.linspace(0,1,num=beta_bins+1,endpoint=True)
ntot=np.linspace(0,max_nbab_to_plot,num=nbab_bins+1,endpoint=True)
e_over_n=np.linspace(0,max_e_over_n,num=e_over_n_bins+1,endpoint=True)
sT3=np.linspace(0,max_sT3,num=sT3_bins+1,endpoint=True)
tt_ov=np.linspace(0,tsize*ovf,num=tcells*ovf+1,endpoint=True)
xx_ov=np.linspace(-xsize*ovf,xsize*ovf,num=xcells*ovf+1,endpoint=True)
yy_ov=np.linspace(-ysize*ovf,ysize*ovf,num=ycells*ovf+1,endpoint=True)
zz_ov=np.linspace(-zsize*ovf,zsize*ovf,num=zcells*ovf+1,endpoint=True)


ptarr=np.linspace(0,ptmax,num=ptbins+1,endpoint=True)
dpt=ptmax/ptbins
raparr=np.linspace(-rapmax,rapmax,num=rapbins+1,endpoint=True)
drap=2*rapmax/rapbins
dbeta=1./beta_bins

part_indx=np.zeros(ntype,dtype=int)
part_indx_raplim=np.zeros(ntype,dtype=int)


#now we fill the arrays with data selected on the particle type
for i in range(good):
    pindex=int(partdata[i,0])
    if(pindex>=0):
      partdata_det[pindex][part_indx[pindex],:]=partdata[i,1:]
      part_indx[pindex]=part_indx[pindex]+1
      if(abs(partdata[i,kbrap])<raplim):
        partdata_det_raplim[pindex][part_indx_raplim[pindex],:]=partdata[i,1:]
        part_indx_raplim[pindex]=part_indx_raplim[pindex]+1

thist_cp=[]
xhist_cp=[]
yhist_cp=[]
zhist_cp=[]
thist=[]
xhist=[]
yhist=[]
zhist=[]
temphist=[]
muhist=[]
betahist=[]
ntothist=[]
e_over_nhist=[]
sT3hist=[]
phase=[]
Tmupoints=[]
Trap=[]
murap=[]
muoverTrap=[]
Tpt=[]
mupt=[]
muoverTpt=[]

print("Now, computing plot data")

for i in range(ntype):
    thist_cp.append(np.histogram(partdata_det[i][:,kt], bins=tt))
    xhist_cp.append(np.histogram(partdata_det[i][:,kx], bins=xx))
    yhist_cp.append(np.histogram(partdata_det[i][:,ky], bins=yy))
    zhist_cp.append(np.histogram(partdata_det[i][:,kz], bins=zz))
    thist.append(np.histogram(partdata_det_raplim[i][:,kt], bins=tt))
    xhist.append(np.histogram(partdata_det_raplim[i][:,kx], bins=xx))
    yhist.append(np.histogram(partdata_det_raplim[i][:,ky], bins=yy))
    zhist.append(np.histogram(partdata_det_raplim[i][:,kz], bins=zz))
    temphist.append(np.histogram(partdata_det_raplim[i][:,ktemp], bins=tb))
    muhist.append(np.histogram(partdata_det_raplim[i][:,kmu], bins=mb))
    betahist.append(np.histogram(partdata_det_raplim[i][:,kbeta], bins=bb))
    ntothist.append(np.histogram(partdata_det_raplim[i][:,knbab], bins=ntot))
    e_over_nhist.append(np.histogram(partdata_det_raplim[i][:,keovern], bins=e_over_n))
    sT3hist.append(np.histogram(partdata_det_raplim[i][:,ksT3], bins=sT3))
    phase.append(np.histogram2d(partdata_det_raplim[i][:,kmu],partdata_det_raplim[i][:,ktemp],bins=(mb,tb)))
    T_avg_p1=0.
    T_err_p1=0.
    mu_avg_p1=0.
    mu_err_p1=0.
    beta_avg_p1=0.
    beta_err_p1=0.
    elements_p1=0.
    T_std_sum_p1=0.
    mu_std_sum_p1=0.
    beta_std_sum_p1=0.
    T_avg_p23=np.zeros(rapbins)
    T_err_p23=np.zeros(rapbins)
    muoverT_avg_p23=np.zeros(rapbins)
    muoverT_err_p23=np.zeros(rapbins)
    mu_avg_p23=np.zeros(rapbins)
    mu_err_p23=np.zeros(rapbins)
    elements_p23=np.zeros(rapbins)
    T_std_sum_p23=np.zeros(rapbins)
    muoverT_std_sum_p23=np.zeros(rapbins)
    mu_std_sum_p23=np.zeros(rapbins)
    T_avg_p45=np.zeros(ptbins)
    T_err_p45=np.zeros(ptbins)
    muoverT_avg_p45=np.zeros(ptbins)
    muoverT_err_p45=np.zeros(ptbins)
    mu_avg_p45=np.zeros(ptbins)
    mu_err_p45=np.zeros(ptbins)
    elements_p45=np.zeros(ptbins)
    T_std_sum_p45=np.zeros(ptbins)
    mu_std_sum_p45=np.zeros(ptbins)
    muoverT_std_sum_p45=np.zeros(ptbins)
    elements_beta=np.zeros(beta_bins)
    if(pcount[i]>1):
      for n in range(pcount[i]):
        if(abs(partdata_det[i][n,krap])<raplim):
          T_avg_p1=T_avg_p1+partdata_det[i][n,ktemp]
          mu_avg_p1=mu_avg_p1+partdata_det[i][n,kmu]
          vbeta=math.sqrt(partdata_det[i][n,kvx]**2+partdata_det[i][n,kvy]**2)
          beta_avg_p1=beta_avg_p1+vbeta
          elements_p1=elements_p1+1
          indx_pt=int(math.floor(partdata_det[i][n,kpt]/dpt))
          T_avg_p45[indx_pt]=T_avg_p45[indx_pt]+partdata_det[i][n,ktemp]
          muoverT_avg_p45[indx_pt]=muoverT_avg_p45[indx_pt]+partdata_det[i][n,ktemp]
          mu_avg_p45[indx_pt]=mu_avg_p45[indx_pt]+partdata_det[i][n,kmu]/partdata_det[i][n,ktemp]
          elements_p45[indx_pt]=elements_p45[indx_pt]+1
        indx_rap=int(math.floor((partdata_det[i][n,krap]+rapmax)/drap))
        T_avg_p23[indx_rap]=T_avg_p23[indx_rap]+partdata_det[i][n,ktemp]
        mu_avg_p23[indx_rap]=mu_avg_p23[indx_rap]+partdata_det[i][n,kmu]
        muoverT_avg_p23[indx_rap]=muoverT_avg_p23[indx_rap]+partdata_det[i][n,kmu]/partdata_det[i][n,kmu]
        elements_p23[indx_rap]=elements_p23[indx_rap]+1
      if(elements_p1 > 1):
        T_avg_p1=T_avg_p1/elements_p1
        mu_avg_p1=mu_avg_p1/elements_p1
        beta_avg_p1=beta_avg_p1/elements_p1
      for h in range(ptbins):
          if(elements_p45[h]>1):
            T_avg_p45[h]=T_avg_p45[h]/elements_p45[h]
            mu_avg_p45[h]=mu_avg_p45[h]/elements_p45[h]
            muoverT_avg_p45[h]=muoverT_avg_p45[h]/elements_p45[h]
      for h in range(rapbins):
          if(elements_p23[h]>1):
            T_avg_p23[h]=T_avg_p23[h]/elements_p23[h]
            mu_avg_p23[h]=mu_avg_p23[h]/elements_p23[h]
            muoverT_avg_p23[h]=muoverT_avg_p23[h]/elements_p23[h]
      for n in range(pcount[i]):
        if(abs(partdata_det[i][n,krap])<raplim):
          T_std_sum_p1=T_std_sum_p1+(partdata_det[i][n,ktemp]-T_avg_p1)**2
          mu_std_sum_p1=mu_std_sum_p1+(partdata_det[i][n,kmu]-mu_avg_p1)**2
          vbeta=math.sqrt(partdata_det[i][n,kvx]**2+partdata_det[i][n,kvy]**2)
          beta_std_sum_p1=beta_std_sum_p1+(vbeta-beta_avg_p1)**2
          indx_pt=int(math.floor(partdata_det[i][n,kpt]/dpt))
          T_std_sum_p45[indx_pt]=T_std_sum_p45[indx_pt]+(partdata_det[i][n,ktemp]-T_avg_p45[indx_pt])**2
          mu_std_sum_p45[indx_pt]=mu_std_sum_p45[indx_pt]+(partdata_det[i][n,kmu]-mu_avg_p45[indx_pt])**2
          muoverT_std_sum_p45[indx_pt]=muoverT_std_sum_p45[indx_pt]+(partdata_det[i][n,kmu]/partdata_det[i][n,ktemp]-muoverT_avg_p45[indx_pt])**2
        indx_rap=int(math.floor((partdata_det[i][n,krap]+rapmax)/drap))
        T_std_sum_p23[indx_rap]=T_std_sum_p23[indx_rap]+(partdata_det[i][n,ktemp]-T_avg_p23[indx_rap])**2
        mu_std_sum_p23[indx_rap]=mu_std_sum_p23[indx_rap]+(partdata_det[i][n,kmu]-mu_avg_p23[indx_rap])**2
        muoverT_std_sum_p23[indx_rap]=muoverT_std_sum_p23[indx_rap]+(partdata_det[i][n,kmu]/partdata_det[i][n,ktemp]-muoverT_avg_p23[indx_rap])**2

      if(elements_p1>1):
        T_err_p1=math.sqrt(T_std_sum_p1/(elements_p1-1))
        mu_err_p1=math.sqrt(mu_std_sum_p1/(elements_p1-1))
        beta_err_p1=math.sqrt(beta_std_sum_p1/(elements_p1-1))
      for h in range(ptbins):
          if(elements_p45[h]>1):
            T_err_p45[h]=math.sqrt(T_std_sum_p45[h]/(elements_p45[h]-1))
            mu_err_p45[h]=math.sqrt(mu_std_sum_p45[h]/(elements_p45[h]-1))
            muoverT_err_p45[h]=math.sqrt(muoverT_std_sum_p45[h]/(elements_p45[h]-1))
      for h in range(rapbins):
          if(elements_p23[h]>1):
            T_err_p23[h]=math.sqrt(T_std_sum_p23[h]/(elements_p23[h]-1))
            mu_err_p23[h]=math.sqrt(mu_std_sum_p23[h]/(elements_p23[h]-1))
            muoverT_err_p23[h]=math.sqrt(muoverT_std_sum_p23[h]/(elements_p23[h]-1))
       
    Tmupoints.append((T_avg_p1,T_err_p1,mu_avg_p1,mu_err_p1,beta_avg_p1,beta_err_p1,elements_p1))
    Tpt.append((T_avg_p45,T_err_p45,elements_p45))
    mupt.append((mu_avg_p45,mu_err_p45,elements_p45))
    muoverTpt.append((muoverT_avg_p45,muoverT_err_p45,elements_p45))
    Trap.append((T_avg_p23,T_err_p23,elements_p23)) 
    murap.append((mu_avg_p23,mu_err_p23,elements_p23)) 
    muoverTrap.append((muoverT_avg_p23,muoverT_err_p23,elements_p23)) 

#now we examine all the particles together
thist_cp.append(np.histogram(partdata[:,kbt], bins=tt))
xhist_cp.append(np.histogram(partdata[:,kbx], bins=xx))
yhist_cp.append(np.histogram(partdata[:,kby], bins=yy))
zhist_cp.append(np.histogram(partdata[:,kbz], bins=zz))
thist.append(np.histogram(partdata_raplim[:,kbt], bins=tt))
xhist.append(np.histogram(partdata_raplim[:,kbx], bins=xx))
yhist.append(np.histogram(partdata_raplim[:,kby], bins=yy))
zhist.append(np.histogram(partdata_raplim[:,kbz], bins=zz))
thist_big=np.histogram(partdata_big[:,kbt], bins=tt)
xhist_big=np.histogram(partdata_big[:,kbx], bins=xx)
yhist_big=np.histogram(partdata_big[:,kby], bins=yy)
zhist_big=np.histogram(partdata_big[:,kbz], bins=zz)
thist_big_ov=np.histogram(partdata_big[:,kbt], bins=tt_ov)
xhist_big_ov=np.histogram(partdata_big[:,kbx], bins=xx_ov)
yhist_big_ov=np.histogram(partdata_big[:,kby], bins=yy_ov)
zhist_big_ov=np.histogram(partdata_big[:,kbz], bins=zz_ov)
temphist.append(np.histogram(partdata_raplim[:,kbtemp], bins=tb))
muhist.append(np.histogram(partdata_raplim[:,kbmu], bins=mb))
betahist.append(np.histogram(partdata_raplim[:,kbbeta], bins=bb))
ntothist.append(np.histogram(partdata_raplim[:,kbnbab], bins=ntot))
e_over_nhist.append(np.histogram(partdata_raplim[:,kbeovern], bins=e_over_n))
sT3hist.append(np.histogram(partdata_raplim[:,kbsT3], bins=sT3))
phase.append(np.histogram2d(partdata_raplim[:,kbmu],partdata_raplim[:,kbtemp],bins=(mb,tb)))
T_avg_p1=0.
T_err_p1=0.
mu_avg_p1=0.
mu_err_p1=0.
beta_avg_p1=0.
beta_err_p1=0.
elements_p1=0.
T_std_sum_p1=0.
mu_std_sum_p1=0.
beta_std_sum_p1=0.
T_avg_p23=np.zeros(rapbins)
T_err_p23=np.zeros(rapbins)
mu_avg_p23=np.zeros(rapbins)
mu_err_p23=np.zeros(rapbins)
muoverT_avg_p23=np.zeros(rapbins)
muoverT_err_p23=np.zeros(rapbins)
elements_p23=np.zeros(rapbins)
T_std_sum_p23=np.zeros(rapbins)
mu_std_sum_p23=np.zeros(rapbins)
muoverT_std_sum_p23=np.zeros(rapbins)
T_avg_p45=np.zeros(ptbins)
T_err_p45=np.zeros(ptbins)
mu_avg_p45=np.zeros(ptbins)
mu_err_p45=np.zeros(ptbins)
muoverT_avg_p45=np.zeros(ptbins)
muoverT_err_p45=np.zeros(ptbins)
elements_p45=np.zeros(ptbins)
T_std_sum_p45=np.zeros(ptbins)
mu_std_sum_p45=np.zeros(ptbins)
muoverT_std_sum_p45=np.zeros(ptbins)
for n in range(good):
  if(abs(partdata[n,kbrap])<raplim):
      T_avg_p1=T_avg_p1+partdata[n,kbtemp]
      mu_avg_p1=mu_avg_p1+partdata[n,kbmu]
      vbeta=math.sqrt(partdata[n,kbvx]**2+partdata[n,kbvy]**2)
      beta_avg_p1=beta_avg_p1+vbeta
      elements_p1=elements_p1+1
      indx_pt=int(math.floor(partdata[n,kbpt]/dpt))
      T_avg_p45[indx_pt]=T_avg_p45[indx_pt]+partdata[n,kbtemp]
      mu_avg_p45[indx_pt]=mu_avg_p45[indx_pt]+partdata[n,kbmu]
      muoverT_avg_p45[indx_pt]=muoverT_avg_p45[indx_pt]+partdata[n,kbmu]/partdata[n,kbtemp]
      elements_p45[indx_pt]=elements_p45[indx_pt]+1
  indx_rap=int(math.floor((partdata[n,kbrap]+rapmax)/drap))
  T_avg_p23[indx_rap]=T_avg_p23[indx_rap]+partdata[n,kbtemp]
  mu_avg_p23[indx_rap]=mu_avg_p23[indx_rap]+partdata[n,kbmu]
  muoverT_avg_p23[indx_rap]=muoverT_avg_p23[indx_rap]+partdata[n,kbmu]/partdata[n,kbtemp]
  elements_p23[indx_rap]=elements_p23[indx_rap]+1

if(elements_p1>1):
  T_avg_p1=T_avg_p1/elements_p1
  mu_avg_p1=mu_avg_p1/elements_p1
  beta_avg_p1=beta_avg_p1/elements_p1
for h in range(ptbins):
    if(elements_p45[h]>1):
        T_avg_p45[h]=T_avg_p45[h]/elements_p45[h]
        mu_avg_p45[h]=mu_avg_p45[h]/elements_p45[h]
        muoverT_avg_p45[h]=muoverT_avg_p45[h]/elements_p45[h]
for h in range(rapbins):
    if(elements_p23[h]>1):
        T_avg_p23[h]=T_avg_p23[h]/elements_p23[h]
        mu_avg_p23[h]=mu_avg_p23[h]/elements_p23[h]
        muoverT_avg_p23[h]=muoverT_avg_p23[h]/elements_p23[h]
for n in range(good):
 if((abs(partdata[n,kbrap])<rapmax) and (partdata[n,kbpt]<ptmax)):
   if(abs(partdata[n,kbrap])<raplim):
       T_std_sum_p1=T_std_sum_p1+(partdata[n,kbtemp]-T_avg_p1)**2
       mu_std_sum_p1=mu_std_sum_p1+(partdata[n,kbmu]-mu_avg_p1)**2
       vbeta=math.sqrt(partdata[n,kbvx]**2+partdata[n,kbvy]**2)
       beta_std_sum_p1=beta_std_sum_p1+(vbeta-beta_avg_p1)**2
       indx_pt=int(math.floor(partdata[n,kbpt]/dpt))
       T_std_sum_p45[indx_pt]=T_std_sum_p45[indx_pt]+(partdata[n,kbtemp]-T_avg_p45[indx_pt])**2
       mu_std_sum_p45[indx_pt]=mu_std_sum_p45[indx_pt]+(partdata[n,kbmu]-mu_avg_p45[indx_pt])**2
       muoverT_std_sum_p45[indx_pt]=muoverT_std_sum_p45[indx_pt]+(partdata[n,kbmu]/partdata[n,kbtemp]-muoverT_avg_p45[indx_pt])**2
   indx_rap=int(math.floor((partdata[n,kbrap]+rapmax)/drap))
   T_std_sum_p23[indx_rap]=T_std_sum_p23[indx_rap]+(partdata[n,kbtemp]-T_avg_p23[indx_rap])**2
   mu_std_sum_p23[indx_rap]=mu_std_sum_p23[indx_rap]+(partdata[n,kbmu]-mu_avg_p23[indx_rap])**2
   muoverT_std_sum_p23[indx_rap]=muoverT_std_sum_p23[indx_rap]+(partdata[n,kbmu]/partdata[n,kbtemp]-muoverT_avg_p23[indx_rap])**2

if(elements_p1>1):
  T_err_p1=math.sqrt(T_std_sum_p1/(elements_p1-1))
  mu_err_p1=math.sqrt(mu_std_sum_p1/(elements_p1-1))
  beta_err_p1=math.sqrt(beta_std_sum_p1/(elements_p1-1))
for h in range(ptbins):
    if(elements_p45[h]>1):
        T_err_p45[h]=math.sqrt(T_std_sum_p45[h]/(elements_p45[h]-1))
        mu_err_p45[h]=math.sqrt(mu_std_sum_p45[h]/(elements_p45[h]-1))
        muoverT_err_p45[h]=math.sqrt(muoverT_std_sum_p45[h]/(elements_p45[h]-1))
for h in range(rapbins):
    if(elements_p23[h]>1):
        T_err_p23[h]=math.sqrt(T_std_sum_p23[h]/(elements_p23[h]-1))
        mu_err_p23[h]=math.sqrt(mu_std_sum_p23[h]/(elements_p23[h]-1))
        muoverT_err_p23[h]=math.sqrt(muoverT_std_sum_p23[h]/(elements_p23[h]-1))
       
Tmupoints.append((T_avg_p1,T_err_p1,mu_avg_p1,mu_err_p1,beta_avg_p1,beta_err_p1,elements_p1))
Tpt.append((T_avg_p45,T_err_p45,elements_p45))
mupt.append((mu_avg_p45,mu_err_p45,elements_p45))
muoverTpt.append((muoverT_avg_p45,muoverT_err_p45,elements_p45))
Trap.append((T_avg_p23,T_err_p23,elements_p23)) 
murap.append((mu_avg_p23,mu_err_p23,elements_p23)) 
muoverTrap.append((muoverT_avg_p23,muoverT_err_p23,elements_p23)) 


print("Pickling (ptype,pcount,pcount_raplim,raplim,ptlim,mb,tb,ptarr,raparr,thist_big,xhist_big,yhist_big,zhist_big,thist_big_ov,xhist_big_ov,yhist_big_ov,zhist_big_ov,thist_cp,xhist_cp,yhist_cp,zhist_cp,thist,xhist,yhist,zhist,temphist,muhist,phase,Tmupoints,Tpt,mupt,muoverTpt,Trap,murap,muoverTrap) as tuples of histograms")
with open(outputfile+"_plot_data_binary.pickle","wb") as po:
     pickle.dump((ptype,pcount,pcount_raplim,raplim,ptlim,mb,tb,ptarr,raparr,thist_big,xhist_big,yhist_big,zhist_big,thist_big_ov,xhist_big_ov,yhist_big_ov,zhist_big_ov,thist_cp,xhist_cp,yhist_cp,zhist_cp,thist,xhist,yhist,zhist,temphist,muhist,betahist,phase,Tmupoints,Tpt,mupt,muoverTpt,Trap,murap,muoverTrap,ntothist,e_over_nhist,sT3hist),po)
print("Printing text data for plots with results from different simulations")
sp="    "
po=open(outputfile+"_point_data.txt","w")
po.write("#Particle    T   T_error    mu     mu_error     <beta>     <beta>_err\n")
for i,v in ptype.items():
    indx=v[0]
    po.write(v[1]+sp+'{:7.3f}'.format(Tmupoints[indx][0])+sp+'{:7.3f}'.format(Tmupoints[indx][1])+sp+'{:7.3f}'.format(Tmupoints[indx][2])+sp+'{:7.3f}'.format(Tmupoints[indx][3])+sp+'{:7.3f}'.format(Tmupoints[indx][4])+sp+'{:7.3f}'.format(Tmupoints[indx][5])+"\n")
po.write("All"+sp+'{:7.3f}'.format(Tmupoints[-1][0])+sp+'{:7.3f}'.format(Tmupoints[-1][1])+sp+'{:7.3f}'.format(Tmupoints[-1][2])+sp+'{:7.3f}'.format(Tmupoints[-1][3])+sp+'{:7.3f}'.format(Tmupoints[-1][4])+sp+'{:7.3f}'.format(Tmupoints[-1][5])+"\n")
po.close()
print("All done.")