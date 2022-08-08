# Author: Gabriele Inghirami - gabriele.inghirami@gmail.com - License GPL v3
import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip

hc3_conv=197.326**3 

#if(len(sys.argv)!=3):
#   print ('Syntax: ./quick_avg.py <inputfile> <energy density limit [GeV/fm^3]>')
#   sys.exit(1)
if(len(sys.argv)!=3):
   print ('Syntax: ./quick_avg.py <inputfile> <min temperature [MeV]>')
   sys.exit(1)

inputfile=sys.argv[1]
lim_temp=np.float64(sys.argv[2])
#limit=np.float64(sys.argv[2])

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  infile=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  infile=open(inputfile,"r")


#count_lines=0
hadrons=0
discarded=0
mu_sum=np.float64(0.)
temp_sum=np.float64(0.)
s_sum=np.float64(0.)
en_sum=np.float64(0.)
nb_sum=np.float64(0.)
tb_sum=np.float64(0.)
part_sum=np.float64(0.)
sT3_ratio_sum=np.float64(0.)
tlow_sum=np.float64(0.)
for line in infile:
#    count_lines=count_lines+1
    stuff=line.split()
#    if(np.float64(stuff[12])>limit):
#      discarded=discarded+1
#      continue
    hadrons=hadrons+1
    temperature=np.float64(stuff[10])
    sT3_ratio=np.float64(stuff[14])*hc3_conv/(temperature**3)
    if(temperature < lim_temp):
      discarded=discarded+1
    else:
      mu_sum=mu_sum+np.float64(stuff[9])
      temp_sum=temp_sum+np.float64(stuff[10])
      nb_sum=nb_sum+np.float64(stuff[11])
      tb_sum=tb_sum+np.float64(stuff[15])
      en_sum=en_sum+np.float64(stuff[12])
      part_sum=part_sum+np.float64(stuff[13])
      s_sum=s_sum+np.float64(stuff[14])
      sT3_ratio_sum=sT3_ratio_sum+sT3_ratio


infile.close()
#total=discarded+hadrons
#print("Energy density limit: "+'{:7.3f}'.format(limit))
#print("Lines read: "+'{:7.3f}'.format(count_lines))
#print("Hadrons: "+str(hadrons)+"  ("+'{:7.3f}'.format(hadrons/total*100)+"%)")
print("Hadrons: "+str(hadrons))
#print("Discarded: "+str(discarded)+"  ("+'{:7.3f}'.format(discarded/total*100)+"%)")
good=hadrons-discarded
print("Accepted: "+str(good)+"  ("+'{:7.3f}'.format(good/hadrons*100)+"%)")
print("Discarded: "+str(discarded)+"  ("+'{:7.3f}'.format(discarded/hadrons*100)+"%)")
print("Average temperature (MeV):  "+'{:7.3f}'.format(temp_sum/good))
print("Average baryon chem. pot. (MeV):  "+'{:7.3f}'.format(mu_sum/good))
print("Average energy density: (GeV/fm^3):  "+'{:7.3f}'.format(en_sum/good))
print("Average net baryon density: (1/fm^3):  "+'{:7.3f}'.format(nb_sum/good))
print("Average total baryon density: (1/fm^3):  "+'{:7.3f}'.format(tb_sum/good))
print("Average particle density (1/fm^3):  "+'{:7.3f}'.format(part_sum/good))
print("Average energy/particle (GeV):  "+'{:7.3f}'.format(en_sum/part_sum))
print("Average entropy density (1/fm^3):  "+'{:7.3f}'.format(s_sum/good))
print("Average s/T^3:  "+'{:7.3f}'.format(sT3_ratio_sum/good))
