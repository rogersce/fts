#!/usr/bin/env python

from numpy import *
import pylo
import visit
import sys
import os
from glob import glob

filelist = sys.argv[1:]
zipped = False

name,sep,ext = filelist[0].rpartition('.')
if ext == 'lzma':
   os.system('tar -xaf '+filelist[0])
   filelist = glob('*.silo')
   zipped = True

#create a .visit file
filelist.sort()
filelist = [fname+'\n' for fname in filelist]
fid = open('db.visit','w')
fid.writelines(filelist)
fid.close()

#plot the database
visit.Launch()
visit.OpenDatabase('db.visit')
visit.AddPlot("Pseudocolor","/psmFpeSolver/data/phi")
visit.DrawPlots()
visit.OpenGUI()

raw_input("Press any key to exit...")

#cleanup
os.system('rm db.visit')
if zipped: os.system('rm *.silo')
