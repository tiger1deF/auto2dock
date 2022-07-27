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
    print(replicateNum)
    
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
GSDock = configBools[10]

if CSVprepare == False:
    if GSDock == True:
        os.system('python3 DockGS.py')
    else:
        os.system('python3 Dock.py')
    sys.exit()


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
AIBindpath = '/home/' + user + '/sars-busters-consolidated/GitCode'
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

def submit_id_mapping(fromDB, toDB, ids):
    API_URL = "https://rest.uniprot.org"
    r = requests.post(
        f"{API_URL}/idmapping/run", data={"from": fromDB, "to": toDB, "ids": ids},
    )
    r.raise_for_status()
    return r.json()["jobId"]


def get_id_mapping_results(job_id):
    POLLING_INTERVAL = 1
    API_URL = "https://rest.uniprot.org"
    while True:
        r = requests.get(f"{API_URL}/idmapping/status/{job_id}")
        r.raise_for_status()
        job = r.json()
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(job["jobStatus"])
        else:
            with open('temp.txt', 'w') as f:
                f.write(json.dumps(job, indent=2))
            with open('temp.txt', 'r') as f:
                lines = f.readlines()
                
            protCounter = 0
            firstConvert = False
            for line in lines:
                if "to" in line:
                    protCounter = protCounter + 1
                    if protCounter == 1:
                        protein =  ast.literal_eval(line[12:])
                    else:
                        if firstConvert == False:
                            tempVal = ast.literal_eval(line[12:])
                            protein = [protein]
                            firstConvert = True
                            protein.append(tempVal)
                        else:
                            protein.append(ast.literal_eval(line[12:]))
            os.remove('temp.txt')
            return protein


def fillarrays():  #Fills arrays for creation of .CSV file
    if NDMscreen == False and GSDock == False:
        
        protList = []
        if CSVprepare == True: #Downloads FASTA from uniprot if configuration set to True
            with open(CSVpath, 'w') as c:
                writer = csv.writer(c)
                headline = ['Symbol','InChiKey', 'SMILE', 'Sequence', 'ID']
                writer.writerow(headline)

            if docking_range[1] == 'all':
                docking_range[1] = rows - 1

            FASTAlist = [] #list of all FASTA sequences
            for i in range(docking_range[0], docking_range[0] + docking_range[1]):
                alphaflag = False
            #for i in range(100, 110): #for testing purposes
                if absnegativeFormat == False:
                    hgncFound = False
                    hgncIndex = 0
                    while hgncFound == False:
                            
                        if 'hgnc' in sheet.cell_value(0, hgncIndex):
                            protnames[i] = sheet.cell_value(i + 1, hgncIndex)
                            
                            hgncFound = True
                        hgncIndex = hgncIndex + 1

                    smiFound = False
                    smiIndex = 0
                    while smiFound == False:
                        if 'SMILE' in sheet.cell_value(0, smiIndex):
                            ligSMILES[i] = sheet.cell_value(i + 1, smiIndex)  #3 for absolute negatives, 4 for AIBIND format
                            smiFound = True
                        smiIndex = smiIndex + 1  #3 for absolute negatives, 4 for AIBIND format
                    try:
                        CIDFound = False
                        CIDindex = 0
                        while CIDFound == False:
                            if 'Chi' in sheet.cell_value(0, CIDindex) or 'Stitch' in sheet.cell_value(0, CIDindex):
                                ligCID[i] = sheet.cell_value(i + 1, CIDindex)
                                CIDFound = True
                                ligID = ligCID[i]
                            CIDindex = CIDindex + 1
                        ligCID[i] = sheet.cell_value(i + 1, ligID)
                    except:
                        ligCID[i] = str(i)
                    if CSVprepare == True: #Gets original inChiKey for AIBind if so configured
                        ligChiKey[i] = ligCID[i]
                    
                    ligCID[i] = (ligCID[i][4:]).lstrip('0') #Processes ligand ID for pcp
                else:
                    ligSMILES[i] = sheet.cell_value(i + 1, 3)
                    ligCID[i] = sheet.cell_value(i + 1, 1)
                    ligChiKey[i] = ligCID[i]
                    colloqNames[i] = sheet.cell_value(i + 1, 4)
                    protnames[i] = colloqNames[i]
                
                #Finds uniprot ID to obtain FASTA sequences
                if absnegativeFormat == True:
                    mapfilename = os.path.join(here, mapfileName)
                    workbook = xlrd.open_workbook(mapfilename)
                    sheets = workbook.sheet_by_index(0)
                    srows = sheets.nrows
                    try:
                        for b in range(0, srows):
                            if colloqNames[i] in str(sheets.cell_value(b + 1, 3)):
                                uniprotID = sheets.cell_value(b + 1, 4)
                                break
                    except:
                        uniprotID = 'N/A'
                    out = mg.querymany(uniprotID, scopes = 'uniprot', fields = 'pdb', species = 'human')
                else:
                    HGNC = protnames[i]  #ID for genecard symbol from excel sheet
                    print(HGNC)
                    out = mg.querymany(HGNC, scopes = 'symbol', fields = 'pdb, uniprot', species = 'human') #Checks the mygene archive
                    #print(out) #for verification purposes
                mygenefailed = False
                downloadfailed = False
                pdbFound = False
                os.chdir(Protpath)
                try:
                    protein = out[0]['pdb']   #Attempts to read PDB ID from uniprot results
                    pdbFound = True
                except KeyError:
                    mygenefailed = True
                    try:
                        uniprot = out[0]['uniprot']['Swiss-Prot']
                        mygenefailed = True
                        alphaFlag = True
                    except:
                        try:
                            os.system('scrapy fetch https://legacy.uniprot.org/uniprot/?query=gene:' + HGNC + ' > out.txt')
                            
                            with open('out.txt', 'r') as f:
                                Data = f.read() 

                            firstResultIndex = Data.index('selectAll-resultStatus')
                            uniprot = Data[firstResultIndex + 58:firstResultIndex + 64]

                            os.remove('out.txt')
                            
                        except Exception as e:
                            print(e)
                            downloadfailed = True

                    try: #remove this block of code if needed
                        job_id = submit_id_mapping(
                            fromDB="UniProtKB_AC-ID", toDB="PDB", ids=[uniprot])
                        protein = get_id_mapping_results(job_id)
                        pdbFound = True
                        mygenefailed = False
                    except:
                        pass

                if mygenefailed == False:
                    if isinstance(protein, list):     #Processes PDB IDs, attempts to download from RCSB
                        RCSBID = [None] * len(protein)
                        for j in range(0, len(RCSBID)):
                            RCSBID[j] = protein[j] + '.pdb'
                    else:
                        RCSBID = protein + '.pdb'
                    try:  #Attempts to download PDB file from RCSB
                        if isinstance(RCSBID, list):
                            c = 0
                            try:
                                for j in range(0, len(RCSBID)):
                                    urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[j], RCSBID[j])
                                    c = c + 1
                                RCSBIDlen = len(RCSBID)
                            except:
                                RCSBIDlen = c
                               
                            size = [None] * RCSBIDlen
                            for j in range(0, RCSBIDlen):
                                with open(RCSBID[j]) as f:
                                    size[j] = len(f.readlines())
                            j = 0
                            d = 0
                            for item in size:   #Filters PDB sequences by size, picks longest sequence
                                if j < item:
                                    j = item
                                    ID = d
                                d = d + 1
                       
                            for j in range(0, RCSBIDlen):
                                if j != ID:
                                    os.remove(RCSBID[j])
                            pfilename = RCSBID[ID]
                            p = pfilename[:-4]
                        else:
                            urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID, RCSBID)
                            pfilename = RCSBID
                            p = pfilename[:-4]
                    except:
                         #Attempts to download and convert mmCIF structure
                        if mmCIF == True and protein == '':
                            
                            if isinstance(protein, list):     #Processes PDB IDs, attempts to download from RCSB
                                
                                RCSBID = [None] * 1 #len(protein)
                                for j in range(0, len(RCSBID)):
                                    RCSBID[j] = protein[j] + '.cif'
                            else:
                                
                                RCSBID = protein + '.cif'
                            try:
                                
                                if isinstance(RCSBID, list):
                                    for j in range(0, len(RCSBID)):
                                        urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[j], RCSBID[j])
                                    size = [None] * len(RCSBID)
                                    for j in range(0, len(RCSBID)):
                                        with open(RCSBID[j]) as f:
                                            size[j] = len(f.readlines())
                                    j = 0
                                    d = 0
                                    for item in size:   #Filters cif sequences by size, picks longest sequence
                                        if j < item:
                                            j = item
                                            ID = d
                                        d = d + 1
                                    for j in range(0, len(RCSBID)):
                                        if j != ID:
                                            os.remove(RCSBID[j])   #Removes all mmCIF files with shorter sequences than the longest file
                                    pfilename = RCSBID[ID]
                                    #urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[0], RCSBID[0])
                                    #pfilename = RCSBID[0]
                                    p = pfilename[:-4]
                                    os.system('gemmi convert ' + pfilename + ' ' + p + '.pdb')   #Converts mmCIF to PDB
                                    os.remove(pfilename)
                                    pfilename = p + '.pdb'
                                    mmCIFflag = True
                                else:
                                
                                    urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID, RCSBID)
                                    pfilename = RCSBID
                                    p = pfilename[:-4]
                                    os.system('gemmi convert ' + pfilename + ' ' + p + '.pdb')   #Converts mmCIF to PDB
                                    os.remove(pfilename)
                                    pfilename = p + '.pdb'
                                    mmCIFflag = True
                            except Exception as e:
                            
                                if absnegativeFormat == False:
                                    out2 = mg.querymany(HGNC, scopes = 'symbol', fields = 'uniprot', species = 'human')
                                    uniprot = out2[0]['uniprot']['Swiss-Prot']
                                else:
                                    uniprot = uniprotID
                               
                                alphaflag = True
                        else:
                            try:
                                #Downloads predicted alphafold structure if PDB ID not recognized
                                if absnegativeFormat == False:
                                    out2 = mg.querymany(HGNC, scopes = 'symbol', fields = 'uniprot', species = 'human')
                                    uniprot = out2[0]['uniprot']['Swiss-Prot']
                                else:
                                    uniprot = uniprotID
                                alphaflag = True
                            except:
                                pass
                else:  #Downloads predicted alphafold structure if no PDB ID found
        
                    try:
                        if absnegativeFormat == False:
                            pfilename = uniprot + '.pdb'
                        else:
                            uniprot = uniprotID

                        alphaflag = True
                    except:
                        print("Protein " + HGNC + " not found in any database!")
                        downloadfailed = True
                    
               

                if alphaflag == True:
                    
                    p = uniprot
                    baseUrl = 'http://www.uniprot.org/uniprot/'
                    currentUrl = baseUrl + uniprot + '.fasta'
                    protList.append(uniprot)
                    
                    response = requests.post(currentUrl)
                    cData = ''.join(response.text)                
                    Seq = StringIO(cData)
                    pSeq = list(SeqIO.parse(Seq,'fasta'))
                    FASTA = [record.seq for record in pSeq]
                    FASTAseq = ''
                    for item in FASTA:
                        FASTAseq = FASTAseq + item
                elif downloadfailed != True:
                    
                    FASTAseq= ''
                    
                    try:
                        out = mg.querymany(HGNC, scopes = 'symbol', fields = 'pdb, uniprot', species = 'human')
                        uniprot = out[0]['uniprot']['Swiss-Prot']
                        baseUrl = 'http://www.uniprot.org/uniprot/'
                        currentUrl = baseUrl + uniprot + '.fasta'
                        protList.append(uniprot)
                        
                        response = requests.post(currentUrl)
                        cData = ''.join(response.text)                
                        Seq = StringIO(cData)
                        pSeq = list(SeqIO.parse(Seq,'fasta'))
                        FASTA = [record.seq for record in pSeq]
                        FASTAseq = ''
                        for item in FASTA:
                            FASTAseq = FASTAseq + item
                        os.remove('out.txt')
                    except:
                        try:
                            os.system('scrapy fetch https://legacy.uniprot.org/uniprot/?query=gene:' + HGNC + ' > out.txt')
                                
                            with open('out.txt', 'r') as f:
                                Data = f.read() 
    
                            firstResultIndex = Data.index('selectAll-resultStatus')
                            uniprot = Data[firstResultIndex + 58:firstResultIndex + 64]
                            baseUrl = 'http://www.uniprot.org/uniprot/'
                            currentUrl = baseUrl + uniprot + '.fasta'
                            protList.append(uniprot)
                            
                            response = requests.post(currentUrl)
                            cData = ''.join(response.text)                
                            Seq = StringIO(cData)
                            pSeq = list(SeqIO.parse(Seq,'fasta'))
                            FASTA = [record.seq for record in pSeq]
                            FASTAseq = ''
                            for item in FASTA:
                                FASTAseq = FASTAseq + item
                            os.remove('out.txt')
                        except:
                            with open(pfilename, 'r') as f:
                                for record in SeqIO.parse(f, 'pdb-atom'):
                                    FASTAseq = FASTAseq + record
                            FASTAseq = FASTAseq.seq          
                          
                else:
                    FASTAseq = '' 
                FASTAlist.append(FASTAseq)

                ID = i 
                with open(CSVpath, 'a') as c:
                    writer = csv.writer(c)
                    if protnames[i] == '' or FASTAseq == '':
                        pass
                    else:
                        row = [protnames[i], ligChiKey[i], ligSMILES[i], FASTAseq, ID]
                        writer.writerow(row)
            os.chdir(Deskpath)
            with open('FASTAList.txt', 'w') as f:
                for item in FASTAlist:
                    f.write(str(item))
    elif GSDock == False:
        pfilename = ''
        uniprot = input("Enter target uniprot ID: ")
        os.chdir(Protpath)
        with open('newprot.txt', 'w') as f:  #For testing purposes
            f.write(uniprot)

        if CSVprepare == True:
            with open(CSVpath, 'w') as c:
                writer = csv.writer(c)
                headline = ['Symbol','InChiKey', 'SMILE', 'Sequence', 'ID']
                writer.writerow(headline)

        
            out = mg.querymany(uniprot, scopes = 'uniprot', fields = 'pdb, uniprot', species = 'human')
            mygenefailed = False
            downloadfailed = False
            pdbFound = False
            os.chdir(Protpath)
            alphaFlag = False
            try:
                protein = out[0]['pdb']   #Attempts to read PDB ID from uniprot results
                pdbFound = True
            except KeyError:                
                try: #remove this block of code if needed
                    job_id = submit_id_mapping(
                        fromDB="UniProtKB_AC-ID", toDB="PDB", ids=[uniprot])
                    protein = get_id_mapping_results(job_id)
                    pdbFound = True
                except Exception as e:
                    pass
                            
            if mygenefailed == False:
                if isinstance(protein, list):     #Processes PDB IDs, attempts to download from RCSB
                    RCSBID = [None] * len(protein)
                    for j in range(0, len(RCSBID)):
                        RCSBID[j] = protein[j] + '.pdb'
                else:
                    RCSBID = protein + '.pdb'
                try:  #Attempts to download PDB file from RCSB
                    if isinstance(RCSBID, list):
                        c = 0
                        try:
                            for j in range(0, len(RCSBID)):
                                urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[j], RCSBID[j])
                                c = c + 1
                            RCSBIDlen = len(RCSBID)
                        except:
                            RCSBIDlen = c
                               
                        size = [None] * RCSBIDlen
                        for j in range(0, RCSBIDlen):
                            with open(RCSBID[j]) as f:
                                size[j] = len(f.readlines())
                        j = 0
                        d = 0
                        
                        for item in size:   #Filters PDB sequences by size, picks longest sequence
                            if j < item:
                                j = item
                                ID = d
                            d = d + 1
                   
                        for j in range(0, RCSBIDlen):
                            if j != ID:
                                os.remove(RCSBID[j])
                        pfilename = RCSBID[ID]
                        p = pfilename[:-4]
                    else:
                        urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID, RCSBID)
                        pfilename = RCSBID
                        p = pfilename[:-4]    
                except:
                     #Attempts to download and convert mmCIF structure
                    if mmCIF == True and pfilename == '':
                        if isinstance(protein, list):     #Processes PDB IDs, attempts to download from RCSB
                            RCSBID = [None] * 1 #len(protein)
                            for j in range(0, len(RCSBID)):
                                RCSBID[j] = protein[j] + '.cif'
                        else:
                            RCSBID = protein + '.cif'
                        try:
                            if isinstance(RCSBID, list):
                                for j in range(0, len(RCSBID)):
                                    urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[j], RCSBID[j])
                                size = [None] * len(RCSBID)
                                for j in range(0, len(RCSBID)):
                                    with open(RCSBID[j]) as f:
                                        size[j] = len(f.readlines())
                                j = 0
                                d = 0
                                for item in size:   #Filters cif sequences by size, picks longest sequence
                                    if j < item:
                                        j = item
                                        ID = d
                                    d = d + 1
                                for j in range(0, len(RCSBID)):
                                    if j != ID:
                                        os.remove(RCSBID[j])   #Removes all mmCIF files with shorter sequences than the longest file
                                pfilename = RCSBID[ID]
                                #urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID[0], RCSBID[0])
                                #pfilename = RCSBID[0]
                                p = pfilename[:-4]
                                os.system('gemmi convert ' + pfilename + ' ' + p + '.pdb')   #Converts mmCIF to PDB
                                os.remove(pfilename)
                                pfilename = p + '.pdb'
                                mmCIFflag = True
                            else:
                                urllib.request.urlretrieve('http://files.rcsb.org/download/' + RCSBID, RCSBID)
                                pfilename = RCSBID
                                p = pfilename[:-4]
                                os.system('gemmi convert ' + pfilename + ' ' + p + '.pdb')   #Converts mmCIF to PDB
                                os.remove(pfilename)
                                pfilename = p + '.pdb'
                                mmCIFflag = True
                        except:
                            if absnegativeFormat == False:
                                out2 = mg.querymany(HGNC, scopes = 'symbol', fields = 'uniprot', species = 'human')
                                uniprot = out2[0]['uniprot']['Swiss-Prot']
                            else:
                                uniprot = uniprotID
                            alphaflag = True
                    else:
                        #Downloads predicted alphafold structure if PDB ID not recognized
                        if absnegativeFormat == False:
                            out2 = mg.querymany(HGNC, scopes = 'symbol', fields = 'uniprot', species = 'human')
                            uniprot = out2[0]['uniprot']['Swiss-Prot']
                        else:
                            uniprot = uniprotID
                           
                        alphaflag = True
            else:  #Downloads predicted alphafold structure if no PDB ID found
                try:
                    if absnegativeFormat == False:
                        pfilename = uniprot + '.pdb'
                    else:
                        uniprot = uniprotID

                    alphaFlag = True
                except:
                    print("Protein " + HGNC + " not found in any database!")
                    downloadfailed = True
                    


            if alphaFlag == True:
                baseUrl = 'http://www.uniprot.org/uniprot/'
                currentUrl = baseUrl + uniprot + '.fasta'
                #protList.append(uniprot)
                response = requests.post(currentUrl)
                cData = ''.join(response.text)                
                Seq = StringIO(cData)
                pSeq = list(SeqIO.parse(Seq,'fasta'))
                FASTA = [record.seq for record in pSeq]
                FASTAseq = ''
                for item in FASTA:
                    FASTAseq = FASTAseq + item
            elif downloadfailed != True:
                FASTAseq= ''
                with open(pfilename, 'r') as f:
                    for record in SeqIO.parse(f, 'pdb-atom'):
                        FASTAseq = FASTAseq + record
                FASTAseq = FASTAseq.seq
                FASTAseq = FASTAseq.replace('X', '')
            else:
                FASTAseq = ''
                
            #calculates HGNC genename for protein
            try:
                job_id = submit_id_mapping(
                    fromDB="UniProtKB_AC-ID", toDB="HGNC", ids=[uniprot])
                protein = get_id_mapping_results(job_id)
                
                HGNC = protein
            except Exception as e:
                try:
                    out = mg.querymany(uniprot, scopes = 'uniprot', fields = 'symbol', species = 'human')
                    HGNC = out[0]['symbol']
                except:
                    print('HGNC ID not found for protein ' + uniprot + '!')
                    HGNC = "N/A"
            for i in range(docking_range[0], endindex):
                smiFound = False
                smiIndex = 0
                while smiFound == False:
                    if 'smi' in sheet.cell_value(0, smiIndex):
                        ligSMILES[i] = sheet.cell_value(i + 1, smiIndex)
                        smiFound = True
                    smiIndex = smiIndex + 1
                    
                keyFound = False
                keyIndex = 0
                while keyFound == False:
                    if 'ick' in sheet.cell_value(0, keyIndex):
                        ligChiKey[i] = sheet.cell_value(i + 1, keyIndex)
                        keyFound = True
                    keyIndex = keyIndex + 1

                ID = i
                with open(CSVpath, 'a') as c:
                    writer = csv.writer(c)
                    if HGNC == '' or FASTAseq == '':
                        row = ['', '', '', '', ID]
                    else:
                        row = [HGNC, ligChiKey[i], ligSMILES[i], FASTAseq, ID]
                    writer.writerow(row) 
            os.chdir(Deskpath)
            with open('FASTAList.txt', 'w') as f:
                for item in FASTAlist:
                    f.write(item)
    else:
        os.chdir(Protpath)
        smiList = []
        protList = []
        ligChiKey = []
        FASTAlist = [] #list of all FASTA sequences
        if CSVprepare == True: #Downloads FASTA from uniprot if configuration set to True
            with open(CSVpath, 'w') as c:
                writer = csv.writer(c)
                headline = ['Symbol','InChiKey', 'SMILE', 'Sequence']
                writer.writerow(headline)
    
            if docking_range[1] == 'all':
                docking_range[1] = rows - 1
            smiFound = False
            smiIndex = 0
            while smiFound == False:
                if 'SMILE' in sheet.cell_value(0, smiIndex) or 'smiles' in sheet.cell_value(0, smiIndex):
                    smiFound = True
                else:
                    smiIndex = smiIndex + 1  #3 for absolute negatives, 4 for AIBIND format
                       
            hgncFound = False
            hgncIndex = 0
            while hgncFound == False:
                if 'hgnc' in sheet.cell_value(0, hgncIndex) or 'PDB' in sheet.cell_value(0, hgncIndex):
                    hgncFound = True
                else:
                    hgncIndex = hgncIndex + 1
            
                    
            for i in range(docking_range[0], docking_range[0] + docking_range[1]):
                smiList.append(sheet.cell_value(i + 1, smiIndex))
                protList.append(sheet.cell_value(i + 1, hgncIndex))
                ligChiKey.append(str(i))
    
                #regex = re.compile('[^a-zA-Z]')
            try:
                for j in range(0, len(protList)):
                    try:
                        FASTAseq= ''
                                
                        try:
                            out = mg.querymany(protList[j], scopes = 'pdb', fields = 'uniprot') #, species = 'human')
                            uniprot = out[0]['uniprot']['Swiss-Prot']
                            baseUrl = 'http://www.uniprot.org/uniprot/'
                            currentUrl = baseUrl + uniprot + '.fasta'
                            protList.append(uniprot)
                                    
                            response = requests.post(currentUrl)
                            cData = ''.join(response.text)                
                            Seq = StringIO(cData)
                            pSeq = list(SeqIO.parse(Seq,'fasta'))
                            FASTA = [record.seq for record in pSeq]
                            FASTAseq = ''
                            for item in FASTA:
                                FASTAseq = FASTAseq + item
                        except Exception as e:
                            try:
                                pfilename = protList[j] + '.pdb'
                                urllib.request.urlretrieve('http://files.rcsb.org/download/' + pfilename, pfilename)
                                with open(pfilename, 'r') as f:
                                    for record in SeqIO.parse(f, 'pdb-atom'):
                                        FASTAseq = FASTAseq + record
                                FASTAseq = FASTAseq.seq
                                os.remove(pfilename)
                            except:
                                pass
                    except:
                        pass
                                  
                    FASTAlist.append(FASTAseq)
                        
                    with open(CSVpath, 'a') as c:
                        writer = csv.writer(c)
                        row = [protList[j], ligChiKey[j], smiList[j], FASTAseq]
                        writer.writerow(row)
                    print(j + 1)
                    j = j + 1
            except Exception as e:
                print(e)
                sys.exit()
        os.chdir(Deskpath)
        with open('FASTAList.txt', 'w') as f:
            for item in FASTAlist:
                f.write(str(item))
                

fillarrays()
os.chdir(AIBindpath)
os.system('python3 test.py')
os.chdir(Deskpath)
os.system('python3 DockTest.py')