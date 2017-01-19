
# coding: utf-8
# A simple hack to use python to submit jobs to bsub. It generates .csh scripts that are then sent to the scheduler.
# It creates file lists for all the samples, as well as output directories on eos
# Author: Stepan Obraztsov

import os, re
import commands
import math, time
import sys
import datetime
import subprocess

print 
print 'START'
print 


#Sample list
samples = ['lower','jet80']



queue = "1nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
#interval = 4 # number files to be processed in a single job, take care to split your file so that you run on all files.


#Everything for data samples here
interval = {}
filepaths = {}
filenumbers = {}
macro = {}
trigger = {}
ismc = {}
isvec = {}
dataprefix = {}

interval['lower'] = 10
interval['jet80'] = 10
interval['mb1'] = 10
interval['mb2'] = 10


filenumbers['lower'] = 500
filenumbers['jet80'] = 500
filenumbers['mb1'] = 880
filenumbers['mb2'] = 618

#For testing on small amount of input files
#for z in samples:
#    filenumbers[z] = 50

trigger['jet80'] = '\"jet80\"'
trigger['lower'] = '\"lower\"'
trigger['mb1'] = '\"minbias\"'
trigger['mb2'] = '\"minbias\"'



ismc['jet80'] = 0
ismc['lower'] = 0
ismc['mb1'] = 0
ismc['mb2'] = 0


dataprefix['jet80'] = 'data_'
dataprefix['lower'] = 'data_'
dataprefix['mb1'] = ''
dataprefix['mb2'] = ''


isvec['jet80'] = 1
isvec['lower'] = 1
isvec['mb1'] = 1
isvec['mb2'] = 1



filepaths['lower'] = "root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtLowerJet_AOD_cmssw758p3_May2016_IVFVtx/HighPtLowerJets/crab_pp5Tev_LowerJet_recalibJP_IVFVtx/160507_222117/"
filepaths['jet80'] = "root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/"
filepaths['mb1'] = "root://xrootd.cmsaf.mit.edu//store/user/abaty/transferTargetDirectories/2015pp_MinBias_2/"
filepaths['mb2'] = "root://xrootd.cmsaf.mit.edu//store/user/abaty/transferTargetDirectories/2015pp_MinBias6/"


#Obsolete. Single macro for all type of samples. 
macro['lower'] = 'deriveClosure'
macro['jet80'] = 'deriveClosure'
macro['mb1'] = 'deriveClosure'
macro['mb2'] = 'deriveClosure'


#Create folders on eos with a time stamp folder{Day}{Month}_{Hour}{Minute}
folder_time = datetime.datetime.now().strftime('%d%m_%H%M')
print folder_time
subprocess.call(["/afs/cern.ch/project/eos/installation/cms/bin/eos.select", "mkdir", "/eos/cms/store/user/stepobr/closure"+folder_time+""])
print "folder closure/"+folder_time+" created"
print 'create folders on eos'
for x in samples:
    subprocess.call(["/afs/cern.ch/project/eos/installation/cms/bin/eos.select", "mkdir", "/eos/cms/store/user/stepobr/closure"+folder_time+"/"+x+""])
    print "folder closure/"+folder_time+"/"+x+" created"


os.system("rm -r jobs_closure")
os.system("mkdir jobs_closure")
os.chdir("jobs_closure")

#generate file lists
for x in samples:
    with open('filelist'+x+'.txt', 'w') as f:
        for number in range(1,filenumbers[x]):
			if number < 1000:
				postfix = "0000/HiForestAOD_" + dataprefix[x]
			elif number < 2000 and number > 999:
				postfix = "0001/HiForestAOD_" + dataprefix[x]
			elif number < 3000 and number > 1999:
				postfix = "0002/HiForestAOD_" + dataprefix[x]
			elif number < 4000 and number > 2999:
				postfix = "0003/HiForestAOD_" + dataprefix[x]
			if x == 'mb1':
				postfix = "MinimumBias2_HiForestAOD_"
			if x == 'mb2':
				postfix = "HiForestAOD_"
			f.write(''+filepaths[x]+postfix+ str(number)+'.root\n')
# generate .csh scripts
for x in samples:
    with open('closure'+x+'.csh', 'w') as fout:
        fout.write("#!/bin/tcsh\n")
        fout.write("setenv HOME /afs/cern.ch/work/s/stepobr/CMSSW_7_5_8/src/HIJetTools/L2Residual\n")
        fout.write("setenv WORK $PWD\n")
        fout.write("setenv X509_USER_PROXY /afs/cern.ch/user/s/stepobr/x509up_u27596\n")
        fout.write("cd ${HOME}\n")
        fout.write("cmsenv\n")
        fout.write("cp "+macro[x]+".C ${WORK}\n")
        fout.write("cp L2Residual.h ${WORK}\n")
        fout.write("cp jobs_closure/filelist"+x+".txt ${WORK}\n")
        fout.write("cp 1801_pythia_response_fits.root ${WORK}\n")
        fout.write("cd ${WORK}\n")
        fout.write("g++ "+macro[x]+".C `root-config --cflags --libs` -O2 -o "+macro[x]+".exe\n")
        fout.write("./"+macro[x]+".exe ${1} ${2} "+str(ismc[x])+" "+str(isvec[x])+" \"filelist"+x+".txt\" "+trigger[x]+"\n")
        fout.write("eos cp test.root /eos/cms/store/user/stepobr/closure"+folder_time+"/"+x+"/ntuple_${1}.root\n")
        fout.write("eos cp "+macro[x]+".C /eos/cms/store/user/stepobr/closure"+folder_time+"/"+macro[x]+".C\n")      
    os.system('chmod 755 closure'+x+'.csh')
   
# submit jobs to bsub
    for y in range(1,filenumbers[x],interval[x]):
    	#Discard the crazy amount of emails sent for every job
        os.system('bsub -o /dev/null -q '+queue+' closure'+x+'.csh '+str(y)+' '+str(y+interval[x]-1)+'')
        #Use the one below for testing
        #os.system('bsub -q '+queue+' closure'+x+'.csh '+str(y)+' '+str(y+interval[x]-1)+'')        
        print "job "+x+" "+str(y)+" submitted"
   
print
#print "your jobs:"
#os.system("bjobs")
print
print 'END'
print



