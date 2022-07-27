#Developed by Christian de Frondeville, with help from Ayan Chatterjee and Michael Sebek
#INSTRUCTIONS

#1. Make sure prepare_receptor4.py, testingfile.py, and Dock.py are in your Linux /home/user directory

#2. Write your Ubuntu username as a string in the user variable definition below
#and the user variable for testing.py (at the top of the code file)

#3. Verify that conda/miniconda is installed, and using pip or conda install, install the following packages:
#pip install xlrd, pip install xlutil, pip install re, pip install bs4, pip install vina, pip install requests
#pip install urllib.request, pip install pubchempy, pip install os, pip install mygene, 
#conda install mgltools, conda install -c conda-forge openbabel, conda install sphinx, conda install numpy
#conda install sphinx_rtd_theme, conda install -c conda-forge boost-cpp, conda install -c conda-forge swig
#conda install gemmi

#4. Verify that the chem-prot-AIBind.xls file is in /home/user directory
#If file is named differently, rename the file in the filename variable for Dock.py and testingfile.py at the top
#If file is formatted differently, edit the columns and rows referenced in testingfile.py

#5. The first parameter in docking_range is the protein you want to start at (where 0 is the first protein)
#and the second parameter is the number of proteins you want to
#Type 'all' for the second parameter to calculate binding energies for all protein-ligand pairs in sheet

#6. Activate your conda env with the downloaded packages 'conda activate <envname>'

#7. Type and enter the command "Python3 Dock.py" in the command line to run the script!

import os
os.environ['MKL_THREADING_LAYER'] = 'GNU'
import sys
import requests
import time
import ast
import json
from urllib.request import urlopen
here = os.path.dirname(os.path.abspath(__file__))  #Prepares Excel file name and path 

destfilename = os.path.join(here, "DockingResults.xls")

with open('dockingRange.txt', 'r') as f:
    dockData = f.readlines()
    if 'all' not in dockData[1][15:]:
        endRange = int(dockData[1][15:])
    else:
        endRange = 'all'
    docking_range = [int(dockData[0][24:]) - 1, endRange] #First parameter is starting index, second is number of proteins indexed
    fname = dockData[2][:-1]
    replicateNum = int(dockData[3][-2:])
    
filename = os.path.join(here, fname)   #Main docking file

#Config options
configBools = []
with open('configuration.txt', 'r') as f:
    configData = f.readlines()
    configLen = len(configData)
user = configData[0][:-1] #Defines program user
for i in range(1, configLen):
    if 'F' in configData[i]:
        configBools.append(False)
    else:
        configBools.append(True)
    
mmCIF = configBools[0]  #Allows code to convert mmCIF files to PDB files        
residueBreak = configBools[3]  #Allows code to break protein up into individual residues
CSVprepare = configBools[6] #Downloads FASTA sequence and adds the sequence and SMILES key into .csv file
absnegativeFormat = configBools[2]
NDMscreen = configBools[8] #Only extracts ligands from respective spreadsheet and dock against a given protein
Alpharemove = configBools[9]
GSDock = configBools[9]

import scrapy   
import shutil 
import re
import statistics
from xlrd import open_workbook
from xlutils.copy import copy  #pip install xlutil
import xlrd
import math
from Bio import SeqIO
wb = xlrd.open_workbook(filename)
sheet = wb.sheet_by_index(0) 
rows = sheet.nrows       #Used by main function as indexer

protnames = [None] * (rows  - 1)#array of protein HGNCs
ligSMILES = [None] * (rows  - 1) #array of ligand identifiers (SMILES)
lignames = [None] * (rows  - 1) #array of ligand names
ligCID = [None] * (rows  - 1) #array of ligand keys
colloqNames = [None] * (rows - 1) #Array of colloquial protein names

mapfileName = 'All_Targets_DrugBank_BindingDB_DTC.xls'
Protpath = '/home/' + user + '/Proteins'  #Paths to all major file locations, user-dependent
Deskpath = '/home/' + user 
ligpath = '/home/' + user + '/Ligands'
indexpath = '/home/' + user + '/index.txt'
vinapath = '/home/' + user + '/vin.txt'
mmcifpath = '/home/' + user + '/mmcif.txt'
outpath = '/home/' + user + '/OutPoses'
predpath = '/home/' + user + '/s4pred'
lookup = 'affinity'
alphafold = 0 #Counts # of alphafold downloads
PDBdown = 0 #Counts # of pdb downloads
faileddownloads = 0 #Counts number of failed downloads

dockFailed = False #tallies if the docking fails in all 3 modes

replicateFailTest = [None] * (rows - 1)  #Prevents program from failing again for additional replicates
if replicateNum > 1:
    dockingenergy = [None] * ((rows - 1) * replicateNum)
else:
    dockingenergy = [None] * (rows - 1)

try:    #Sets up protein, ligand, and output pose folders, creates them if they do not exist
    os.chdir(Protpath)
except:
    os.mkdir(Protpath)
    print('Performing first-time setup of docking environment!')
try:
    os.chdir(ligpath)
except:
    os.mkdir(ligpath)
try:
    os.chdir(outpath)
except:
    os.mkdir(outpath)

if CSVprepare == True: #Initializes needed settings if AIBind used to identify binding locations
    import csv 
    import re
    import mygene  #pip install mygene
    import requests
    from Bio import SeqIO
    from io import StringIO
    from Bio.PDB import PDBList, PDBParser, MMCIFIO #pip install Bio.PDB
    import urllib
    from urllib import request
    import urllib.request    #pip install urllib.request
    mg = mygene.MyGeneInfo()  #instantiates mygene library
    CSVpath = '/home/' + user + '/Proteins/AIBindCSV.csv'
    ligChiKey = [None] * (rows - 1)

if docking_range[1] == 'all':   #Uses docking_range to determine end of runtime
    endindex = rows - 1
else:
    endindex = docking_range[1] + docking_range[0]

def strucPredict():
    strucList = []
    print('Now predicting secondary structure!')
    os.chdir(Deskpath)
    c = 0
    with open('FASTAList.txt', 'r') as f:
        dat = f.readlines()
    FASTAlist = []
    for line in dat:
        FASTAlist.append(line)
    os.remove('FASTAList.txt')
    os.chdir(predpath)
    for item in FASTAlist:
        print(c)
        strucList.append([])
        try:
            with open(str(c) + '.fasta', 'w') as f:
                f.write('>' + str(c) + '\n')
                f.write(str(item))
            os.system('python3 run_model.py ' + str(c)  + '.fasta > ' + str(c) + '1.fasta')
            with open(str(c) + '1.fasta', 'r') as f:
                data = f.readlines()
            os.remove(str(c) + '1.fasta')
            os.remove(str(c) + '.fasta')
        except Exception as e:
            print(e)
            data = ['']
              
        a = 0
        for line in data:
            if a > 1:
                strucList[c].append(line[7:8])
                    #FOR FALSE POSITIVE TESTING
                    #strucList[c].append('C')
            a = a + 1
        c = c + 1
                
    os.chdir(Protpath)
    
    with open('bindingLocations.txt', 'r') as f:
        AIBindsites = ast.literal_eval(f.read())
                
    #AIBindsites = AIBindsites[docking_range[0]:docking_range[0] + docking_range[1]]
    pairIndex = 0
    newSites = []
    breakFlag = False

    for item in AIBindsites:  #Modifies binding sites from AIBind for secondary structure
        newSites.append([])
        newSiteIndex = 0
        for pair in item:
            if pair[1] >= len(strucList[pairIndex]) - 1:
                if pair[0] >= len(strucList[pairIndex]) - 1:
                    pair = [len(strucList[pairIndex]) - 1, len(strucList[pairIndex]) - 1]
                else:
                    pair[1] = len(strucList[pairIndex]) - 1
            newSites[pairIndex].append([])
            breakFlag = False
            if pair[0] == pair[1]:
                pass
            else:
                for c in range(pair[0], pair[1] + 1):
                    if strucList[pairIndex][c] == 'H':
                        if breakFlag == False:
                            breakFlag = True
                            if c != pair[0] and len(newSites[pairIndex][newSiteIndex]) == 1:
                                            
                                newSites[pairIndex][newSiteIndex].append(c)
                                newSites[pairIndex].append([])
                                newSiteIndex = newSiteIndex + 1
                    else:
                        if breakFlag == True:
                                    
                            if c != pair[1] and c != pair[0]:
                                        
                                newSites[pairIndex][newSiteIndex].append(c)
                            breakFlag = False
                        else:
                            if c == pair[0]:
                                #newSites[pairIndex].append([])
                                newSites[pairIndex][newSiteIndex].append(c)
                            elif c == pair[1]:
                                newSites[pairIndex][newSiteIndex].append(c)
                                        
                
            if len(newSites[pairIndex][newSiteIndex]) > 0:
                        newSiteIndex = newSiteIndex + 1
        pairNum = 0
        for pair in newSites[pairIndex]:
            if pair == []:
                del newSites[pairIndex][pairNum]
            else:
                        
                pair[0] = int(pair[0] / 3)
                pair[1] = int(pair[1] / 3)
                pairNum = pairNum + 1
        pairIndex = pairIndex + 1
                                                 
    os.chdir(Protpath)
                    
    with open('bindingSites.txt', 'w') as f: #Eventually write entire list to file 
        f.write(str(newSites))
                
    
strucPredict()

#sys.exit()
def reset():  #empties all directories of generated folders
    if NDMscreen == False:
        os.chdir(Protpath)
        protfiles = os.listdir(Protpath)  #Clears protein folder of protein .pdbqt files
        for item in protfiles:
            if item.endswith(".sdf") or item.endswith(".pdb") or item.endswith(".pdbqt"):
                os.remove(item)

    os.chdir(ligpath) 
    ligfiles = os.listdir(ligpath)   #Clears ligand folder of ligand .pdbqt files
    for item in ligfiles:
        if item.endswith(".pdbqt") or item.endswith(".pdb") or item.endswith(".sdf"):
            os.remove(item)

    os.chdir(Deskpath)
    deskfiles = os.listdir(Deskpath)
    for item in deskfiles:
        if item.endswith(".pdbqt"):
            os.remove(item)

class breakLoop(Exception): pass

def peptideFreqCalc():  #Use to calculate values for peptide binding frequency
    global dockingNums, arrayNum, triNum
    
    if replicateFailTest[i] == None:
        if GSDock == False:
            os.system('python3 TestingFile.py > output.txt') #Writes docking output to file for later reading
        else:
            os.system('python3 testingfileGS.py > output.txt')
        
    with open(mmcifpath, 'r') as f:
        index = int(f.read())

    #Isolate which section of the protein trigrams are calculated, must be changed in testingfile.py 
    #index = 1 #Total number of protein subsections docked
    startindex = 0 #Start of subsection range
    
    with open(indexpath, 'w') as f:  #Writes iteration number to indexer file
        f.write(str(i))
                        
    kcal = []
    dockingNums = []
    modeindex = 0

    try:
        os.chdir(Protpath)
        with open('noPockets.txt', 'r') as f:
            pass
        ans = 'Non-binding' #No binding pockets found
        print('No binding pockets found for ligand!')
    except:
        os.chdir(Deskpath)
        with open('output.txt', 'r') as myFile:
            data = myFile.readlines()
            datalen = len(data)
            os.chdir(Deskpath)
            with open('trigramNum.txt', 'r') as f:
                triNum = int(f.read())
            #os.chdir(Deskpath)
            if triNum == 0:
                raise breakLoop
        
            index = triNum + startindex
            
            for e in range(0, index):
                for b in range(modeindex + 1, datalen):
                    if lookup in data[b]:
                        modeindex = b   #Script starts reading from "affinity" keyword in output.txt results
                        lookedup = True
                        break
                for j in range(0, 8):    #Iterates through all 8 poses
                    try:
                        binding = re.sub('[^0-9^.^-]', '', data[modeindex + 3 + j][5:20])
                        scitest = binding.partition('.')
                        if '-0' in scitest[2]:
                            scitest = list(scitest)
                            exponent = scitest[2].strip(' ')
                            base = float(binding.strip(' ')[:-3])
                            exponent = exponent[-3:]
                            exponent = int(exponent)
                            result = pow(10, exponent) * base
                            kcal.append(abs(float(result)))
                        else:
                            kcal.append(abs(float(binding)))
                    except:
                        True
                dockingNums.append(max(kcal) * -1)
                kcal = []
            ans = min(dockingNums)
    return ans
       
            
    #Create array of numbers in proper range to graph against binding profile
    arrayNum = [i for i in range(startindex + 1, index + 1)]


os.chdir(Deskpath)  #Runs loops multiple times (can average them later)
if docking_range[1] == 'all':
    totalProg = replicateNum * ((rows - 1) - docking_range[0])  #Progress divisor
else:
    totalProg = docking_range[1] * replicateNum
progcount = 0 #Progression indexer
rb = open_workbook(destfilename,formatting_info=True)  #Copies binding affinity to spreadsheet
wb = copy(rb)
sheet = wb.get_sheet(0) 
#with open('dub.txt', 'w') as f: #for testing
    #f.write(str(docking_range[0]) + ' ' + str(endindex) + ' ' + str(replicateNum))
for c in range(0, replicateNum):  #Replicate repeater loop
    print(c)
    with open('output.txt', 'w') as f:
        True   #Clears output file upon moving to next replicate
    for i in range(docking_range[0], endindex):  #Main runtime loop
        os.chdir(Deskpath)
        progress =  str(100 * (progcount / totalProg)) + '%'
        if progress == '100.0%':
            print('Progress: ' + progress)
            break  #Ends runtime loop when progress reaches 100%
        else:
            if i == docking_range[0] and c == 0:
                print('Progress: 0%')   #Prints initial progress value 
            else:
                print('Progress: ' + progress)   #Prints calculated progress value based on replicates and # of pairs
        progcount = progcount + 1
        print(i + 1) #Index of protein/ligand being docked
        
        with open(vinapath, 'w') as f:
            f.write('SUCCESS')
        with open('safeMode.txt', 'w') as f:
            f.write('SUCCESS')
        
        if i > (docking_range[0] - 1) and i < endindex:
            if residueBreak == False:
                with open(indexpath, 'w') as f:  #Writes iteration number to indexer file
                    f.write(str(i))
                if replicateFailTest[i] == None:
                    os.system('python3 testingfileTest.py > output.txt') #Writes docking output to file for later reading
                
                try:    #Processes multiple runs of docking for subsectioned protein
                    with open(mmcifpath, 'r') as f:
                        index = int(f.read())   #Number of subsections read from file
                    kcal = []
                    dockingnums = []
                    modeindex = 0
                    with open('output.txt', 'r') as myFile:
                        data = myFile.readlines()
                    datalen = len(data)
                    for e in range(0, index):
                        for b in range(modeindex + 1, datalen):
                            if lookup in data[b]:
                                modeindex = b   #Script starts reading from "affinity" keyword in output.txt results
                                lookedup = True
                                break
                        for j in range(0, 8):    #Iterates through all 8 poses for mmCIF
                            try:
                                binding = re.sub('[^0-9^.^-]', '', data[modeindex + 3 + j][5:20])
                                scitest = binding.partition('.')
                                if '-0' in scitest[2]:
                                    scitest = list(scitest)
                                    exponent = scitest[2].strip(' ')
                                    base = float(binding.strip(' ')[:-3])
                                    exponent = exponent[-3:]
                                    exponent = int(exponent)
                                    result = pow(10, exponent) * base
                                    kcal.append(abs(float(result)))
                                else:
                                    kcal.append(abs(float(binding)))
                            except:
                                True
                            dockingnums.append(max(kcal))
                            kcal = []
                        max(dockingnums) * -1   
                        dockingenergy[i + ((rows - 1) * c)] = max(dockingnums) * -1
                except:
                    lookedup = False
                    with open('output.txt') as myFile:
                        for num, line in enumerate(myFile, 1):
                            if lookup in line:
                                modeindex = num
                                lookedup = True  #Isolates binding energies from terminal output       
                        if lookedup == False:
                            replicateFailTest[i] = 'FAILED'
                            with open(vinapath, 'w') as f:
                                f.write('FAILED')
                            os.system('python3 testingfileTest.py > output.txt')
                            try:
                                with open(mmcifpath, 'r') as f:
                                    index = int(f.read())   #Number of subsections read from file
                            except:
                                index = 1
                            try:
                                kcal = []
                                dockingnums = []
                                modeindex = 0
                                with open('output.txt', 'r') as myFile:
                                    data = myFile.readlines()
                                datalen = len(data)
                                for e in range(0, index):
                                    for b in range(modeindex + 1, datalen):
                                        if lookup in data[b]:
                                            modeindex = b   #Script starts reading from "affinity" keyword in output.txt results
                                            lookedup = True
                                            break
                                    for j in range(0, 8):    #Iterates through all 8 poses for mmCIF
                                        try:
                                            binding = re.sub('[^0-9^.^-]', '', data[modeindex + j + 3][5:20])
                                            scitest = binding.partition('.')
                                            if '-0' in scitest[2]:
                                                scitest = list(scitest)
                                                exponent = scitest[2].strip(' ')
                                                base = float(binding.strip(' ')[:-3])
                                                exponent = exponent[-3:]
                                                exponent = int(exponent)
                                                result = pow(10, exponent) * base
                                                kcal.append(abs(float(result)))
                                            else:
                                                kcal.append(abs(float(binding)))
                                        except:
                                            True
                                        dockingnums.append(max(kcal))
                                        kcal = []
                            except:
                                print('Running program in safe mode!')
                                with open(vinapath, 'w') as f:
                                    f.write('FAILED')
                                with open('safeMode.txt', 'w') as f:
                                    f.write('FAIL')
                                os.system('python3 testingfileTest.py > output.txt')
                                with open('safeMode.txt', 'w') as f:
                                    f.write('')
                                try:
                                    with open(mmcifpath, 'r') as f:
                                        index = int(f.read())   #Number of subsections read from file
                                except:
                                    index = 1
                                try:
                                    kcal = []
                                    dockingnums = []
                                    modeindex = 0
                                    with open('output.txt', 'r') as myFile:
                                        data = myFile.readlines()
                                    datalen = len(data)
                                    for e in range(0, index):
                                        for b in range(modeindex + 1, datalen):
                                            if lookup in data[b]:
                                                modeindex = b   #Script starts reading from "affinity" keyword in output.txt results
                                                lookedup = True
                                                break
                                        for j in range(0, 8):    #Iterates through all 8 poses for mmCIF
                                            try:
                                                binding = re.sub('[^0-9^.^-]', '', data[modeindex + j + 3][5:20])
                                                scitest = binding.partition('.')
                                                if '-0' in scitest[2]:
                                                    scitest = list(scitest)
                                                    exponent = scitest[2].strip(' ')
                                                    base = float(binding.strip(' ')[:-3])
                                                    exponent = exponent[-3:]
                                                    exponent = int(exponent)
                                                    result = pow(10, exponent) * base
                                                    kcal.append(abs(float(result)))
                                                else:
                                                    kcal.append(abs(float(binding)))
                                            except:
                                                True
                                            dockingnums.append(max(kcal))
                                            kcal = []
                                except:
                                    print('Docking was unable to be performed!')
                                    dockFailed = True
                            if dockFailed != True:
                                dockingenergy[i + ((rows - 1) * c)] = max(dockingnums) * -1
                            else:
                                dockingenergy[i + ((rows - 1) * c)] = 'N/A'
                        else:
                        
                            kcal = []
                            with open('output.txt', 'r') as myFile:
                                data = myFile.readlines()
                                for j in range(0, 20):    #Iterates through all 20 poses
                                    try:
                                        binding = re.sub('[^0-9^.^-]', '', data[modeindex + j + 2][5:20])
                                        scitest = binding.partition('.')
                                        if '-0' in scitest[2]:
                                            scitest = list(scitest)
                                            exponent = scitest[2].strip(' ')
                                            base = float(binding.strip(' ')[:-3])
                                            exponent = exponent[-3:]
                                            exponent = int(exponent)
                                            result = pow(10, exponent) * base
                                            kcal.append(abs(float(result)))
                                        else:
                                            if binding != '':
                                                kcal.append(abs(float(binding)))
                                    except:
                                        True
                                if kcal != []:
                                    dockingenergy[i + ((rows - 1) * c)] = max(kcal) * -1
                                else:
                                    dockingenergy[i + ((rows - 1) * c)] = 'N/A'
                                #print(dockingenergy[i + ((rows - 1) * c)])
            else:
                try:
                    with open(indexpath, 'w') as f:  #Writes iteration number to indexer file
                        f.write(str(i))
                    minDock = peptideFreqCalc()
                    dockingenergy[i + ((rows - 1) * c)] = minDock
                except Exception as e:
                    print(e)
                    try:
                        os.chdir(Protpath)
                        with open('noPockets.txt', 'r') as f:
                            pass
                        dockingenergy[i + ((rows - 1) * c)] = 'Non-binding' #No binding pockets found
                        print('No binding pockets found for ligand!')
                    except Exception as e:
                        print(e)
                        print('Protein Dissection Failed!')
            if NDMscreen == False:
                try:
                    os.chdir(Deskpath)
                    with open('type.txt', 'r') as f:
                        line = f.readlines()
                        if 'ALPHAFOLD' in line[0]:    #testingfile outputs type of download to output.txt, which is read by Dock.py
                            alphafold = alphafold + 1
                        else:
                            PDBdown = PDBdown + 1
                        if dockingenergy[i + ((rows - 1) * c)] == 'N/A':
                            faileddownloads = faileddownloads + 1
                except Exception as e:
                    print(e)
                    pass
                    
        try:
            os.remove(mmcifpath)
        except:
            True
        #reset()
        sheet.write(i + 1, 0 + c, dockingenergy[i + ((rows - 1) * c)])
        wb.save(destfilename)
        try:
            os.chdir(Protpath)
            os.remove('noPockets.txt')
        except:
            pass
        
os.chdir(Deskpath)
try:
    #os.remove('output.txt')  #Text file cleanup
    os.remove('index.txt')
except:
    pass
try:
    os.remove('vin.txt')
except:
    pass
try:
    os.remove('type.txt')
except:
    pass
try:
    os.chdir(Protpath)
    os.remove('protname.txt')
    os.remove('newprot.txt')
except:
    pass

rb = open_workbook(filename,formatting_info=True)  #Copies binding affinity to spreadsheet
wb = copy(rb)
sheet = wb.get_sheet(0)
sheet.write(0, 7, 'Binding Affinity (kcal/mol)')
for j in range(0, replicateNum):
    for i in range(rows - 1):  #writes docking energies back to chemical excel sheet
        sheet.write(i + 1, 7 + j, dockingenergy[i + ((rows - 1) * j)])
if replicateNum > 1:  #Calculates standard deviation for multiple replicates
    faileddownloads = faileddownloads / replicateNum
    alphafold = alphafold / replicateNum
    PDBdown = PDBdown / replicateNum
    sheet.write(0, 7 + replicateNum, 'Standard Deviation')
    for i in range(docking_range[0], docking_range[1]):
        STD = [None] * replicateNum
        for c in range(replicateNum):
            STD[c] = dockingenergy[i + ((rows - 1) * c)]
        sheet.write(i + 1, 7 + replicates, statistics.stdev(STD))
wb.save(filename)  #Saves changes to chem excel sheet

os.chdir(Protpath)
#protfiles = os.listdir(Protpath)  #Clears protein folder of protein .pdbqt files
#for item in protfiles:
    #if item.endswith(".txt") or item.endswith(".pdbqt"):
        #os.remove(item)
                
if NDMscreen == False:
    print('Failed downloads: ' + str(faileddownloads))   
    print('Alphafold estimated configurations downloaded:' + str(alphafold))
    print('Observed crystallized PDB configurations downloaded: ' + str(PDBdown))