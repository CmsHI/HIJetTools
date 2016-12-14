
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
samples = ['lower','jet80','mb1','mb2','pthat30','pthat50','pthat80','pthat120']
#samples = ['lower','jet80','pthat30','pthat50','pthat80','pthat120','herwig50','herwig100','herwig200']
#samples = ['pthat50','pthat80','pthat120','herwig50','herwig100','herwig200']
#samples = ['herwig50','herwig100','herwig200']

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

interval['lower'] = 4
interval['jet80'] = 4
interval['mb1'] = 4
interval['mb2'] = 4
interval['pthat15'] = 4
interval['pthat30'] = 4
interval['pthat50'] = 4
interval['pthat80'] = 4
interval['pthat120'] = 4
interval['herwig50'] = 5
interval['herwig100'] = 5
interval['herwig200'] = 5

filenumbers['lower'] = 500
filenumbers['jet80'] = 500
filenumbers['mb1'] = 880
filenumbers['mb2'] = 618
filenumbers['pthat15'] = 500
filenumbers['pthat30'] = 500
filenumbers['pthat50'] = 500
filenumbers['pthat80'] = 500
filenumbers['pthat120'] = 500
filenumbers['herwig50'] = 1802
filenumbers['herwig100'] = 3251
filenumbers['herwig200'] = 1569

#For testing on small amount of input files
#for z in samples:
#    filenumbers[z] = 10

trigger['jet80'] = '\"jet80\"'
trigger['lower'] = '\"lower\"'
trigger['mb1'] = '\"minbias\"'
trigger['mb2'] = '\"minbias\"'
trigger['pthat15'] = '\"pthat15\"'
trigger['pthat30'] = '\"pthat30\"'
trigger['pthat50'] = '\"pthat50\"'
trigger['pthat80'] = '\"pthat80\"'
trigger['pthat120'] = '\"pthat120\"'
trigger['herwig50'] = '\"herwig50\"'
trigger['herwig100'] = '\"herwig100\"'
trigger['herwig200'] = '\"herwig200\"'


ismc['jet80'] = 0
ismc['lower'] = 0
ismc['mb1'] = 0
ismc['mb2'] = 0
ismc['pthat15'] = 1
ismc['pthat30'] = 1
ismc['pthat50'] = 1
ismc['pthat80'] = 1
ismc['pthat120'] = 1
ismc['herwig50'] = 1
ismc['herwig100'] = 1
ismc['herwig200'] = 1

dataprefix['jet80'] = 'data_'
dataprefix['lower'] = 'data_'
dataprefix['mb1'] = ''
dataprefix['mb2'] = ''
dataprefix['pthat15'] = ''
dataprefix['pthat30'] = ''
dataprefix['pthat50'] = ''
dataprefix['pthat80'] = ''
dataprefix['pthat120'] = ''
dataprefix['herwig50'] = ''
dataprefix['herwig100'] = ''
dataprefix['herwig200'] = ''

isvec['jet80'] = 1
isvec['lower'] = 1
isvec['mb1'] = 1
isvec['mb2'] = 1
isvec['pthat15'] = 0
isvec['pthat30'] = 0
isvec['pthat50'] = 0
isvec['pthat80'] = 0
isvec['pthat120'] = 0
isvec['herwig50'] = 1
isvec['herwig100'] = 1
isvec['herwig200'] = 1


filepaths['lower'] = "root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtLowerJet_AOD_cmssw758p3_May2016_IVFVtx/HighPtLowerJets/crab_pp5Tev_LowerJet_recalibJP_IVFVtx/160507_222117/"
filepaths['jet80'] = "root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/"
filepaths['mb1'] = "root://xrootd.cmsaf.mit.edu//store/user/abaty/transferTargetDirectories/2015pp_MinBias_2/"
filepaths['mb2'] = "root://xrootd.cmsaf.mit.edu//store/user/abaty/transferTargetDirectories/2015pp_MinBias6/"
filepaths['pthat15'] ="root://xrootd.cmsaf.mit.edu//store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet15_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_125350/" 
filepaths['pthat30'] ="root://xrootd.cmsaf.mit.edu//store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133336/"
filepaths['pthat50'] ="root://xrootd.cmsaf.mit.edu//store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet50_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133347/"
filepaths['pthat80'] ="root://xrootd.cmsaf.mit.edu//store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133359/"
filepaths['pthat120'] ="root://xrootd.cmsaf.mit.edu//store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133414/"
filepaths['herwig50'] = "root://xrootd.cmsaf.mit.edu//store/user/dgulhan/QCD_Pt_50to100_TuneEE5C_502TeV_herwigpp_cff_py/HiForest_QCD_Pt50to100_TuneEE5C_502TeV_herwigpp_cff_py/160628_085221/"
filepaths['herwig100'] = "root://xrootd.cmsaf.mit.edu//store/user/dgulhan/QCD_Pt_100to200_TuneEE5C_502TeV_herwigpp_cff_py/HiForest_QCD_Pt100to200_TuneEE5C_502TeV_herwigpp_cff_py/160628_085317/"
filepaths['herwig200'] = "root://xrootd.cmsaf.mit.edu//store/user/dgulhan/QCD_Pt_200_TuneEE5C_502TeV_herwigpp_cff_py/HiForest_QCD_Pt200to9999_TuneEE5C_502TeV_herwigpp_cff_py/160628_090618/"

#Obsolete. Single macro for all type of samples. 
macro['lower'] = 'deriveResponse'
macro['jet80'] = 'deriveResponse'
macro['mb1'] = 'deriveResponse'
macro['mb2'] = 'deriveResponse'
macro['pthat15'] = 'deriveResponse'
macro['pthat30'] = 'deriveResponse'
macro['pthat50'] = 'deriveResponse'
macro['pthat80'] = 'deriveResponse'
macro['pthat120'] = 'deriveResponse'
macro['herwig50'] = 'deriveResponse'
macro['herwig100'] = 'deriveResponse'
macro['herwig200'] = 'deriveResponse'

#Create folders on eos with a time stamp folder{Day}{Month}_{Hour}{Minute}
folder_time = datetime.datetime.now().strftime('%d%m_%H%M')
print folder_time
subprocess.call(["/afs/cern.ch/project/eos/installation/cms/bin/eos.select", "mkdir", "/eos/cms/store/user/stepobr/res"+folder_time+""])
print "folder res/"+folder_time+" created"
print 'create folders on eos'
for x in samples:
    subprocess.call(["/afs/cern.ch/project/eos/installation/cms/bin/eos.select", "mkdir", "/eos/cms/store/user/stepobr/res"+folder_time+"/"+x+""])
    print "folder res/"+folder_time+"/"+x+" created"


os.system("rm -r jobs_response")
os.system("mkdir jobs_response")
os.chdir("jobs_response")

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
    with open('response'+x+'.csh', 'w') as fout:
        fout.write("#!/bin/tcsh\n")
        fout.write("setenv HOME /afs/cern.ch/work/s/stepobr/CMSSW_7_5_8/src/HIJetTools/L2Residual\n")
        fout.write("setenv WORK $PWD\n")
        fout.write("setenv X509_USER_PROXY /afs/cern.ch/user/s/stepobr/x509up_u27596\n")
        fout.write("cd ${HOME}\n")
        fout.write("cmsenv\n")
        fout.write("cp "+macro[x]+".C ${WORK}\n")
        fout.write("cp jobs_response/filelist"+x+".txt ${WORK}\n")
        fout.write("cd ${WORK}\n")
        fout.write("g++ "+macro[x]+".C `root-config --cflags --libs` -O2 -o "+macro[x]+".exe\n")
        fout.write("./"+macro[x]+".exe ${1} ${2} "+str(ismc[x])+" "+str(isvec[x])+" \"filelist"+x+".txt\" "+trigger[x]+"\n")
        fout.write("eos cp test.root /eos/cms/store/user/stepobr/res"+folder_time+"/"+x+"/ntuple_${1}.root\n")
        fout.write("eos cp "+macro[x]+".C /eos/cms/store/user/stepobr/res"+folder_time+"/"+macro[x]+".C\n")      
    os.system('chmod 755 response'+x+'.csh')
   
# submit jobs to bsub
    for y in range(1,filenumbers[x],interval[x]):
    	#Discard the crazy amount of emails sent for every job
        os.system('bsub -o /dev/null -q '+queue+' response'+x+'.csh '+str(y)+' '+str(y+interval[x]-1)+'')
        #Use the one below for testing
        #os.system('bsub -q '+queue+' response'+x+'.csh '+str(y)+' '+str(y+interval[x]-1)+'')        
        print "job "+x+" "+str(y)+" submitted"
   
print
#print "your jobs:"
#os.system("bjobs")
print
print 'END'
print



