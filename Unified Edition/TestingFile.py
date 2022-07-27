import os
import json
here = os.path.dirname(os.path.abspath(__file__))  #Prepares Excel file name and path
with open('dockingRange.txt', 'r') as f:
    dockData = f.readlines()
    if 'all' not in dockData[1][15:]:
        endRange = int(dockData[1][15:])
    else:
        endRange = 'all'
    docking_range = [int(dockData[0][24:]) - 1, endRange] #First parameter is starting index, second is number of proteins indexed
    fname = dockData[2][:-1]  
    
filename = os.path.join(here, fname)

from os import listdir                #NOTE: Requires installation of ADGFRSuite
import ast
import xlrd  #pip install xlrd
import requests as r   #pip install requests
import pubchempy as pcp   #pip install pubchempy
import sys
import urllib
from urllib import request
import urllib.request    #pip install urllib.request
import shutil   #pip install shutil
import openbabel #openbabel 3.1.0 (installed using a .whl script)
from openbabel import pybel
from Bio import PDB #conda/pip install Bio
from Bio import SeqIO
from io import StringIO
from Bio.PDB import PDBList, PDBParser, MMCIFIO #pip install Bio.PDB
import mygene  #pip install mygene
import re
from bs4 import BeautifulSoup #pip install bs4
import statistics
from vina import Vina  #requires sphinx, sphinx_rtd_theme, numpy
import gemmi
import collections
#import pymol2

#Configuration of settings
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
SMILEUse = configBools[1] #Allows code to use SMILES sequences as ligand files
absnegativeFormat = configBools[2] #Allows code to parse data resembling the absnegativeBindingDB excel sheet format
residueBreak = configBools[3]  #Allows code to break protein up into individual residues
trigramCombine = configBools[4] #Allows code to group broken residues into trigrams
HETremove = configBools[5] #If code is broken up into trigrams, removes HETATMS from file before processing
CSVprepare = configBools[6] #Downloads FASTA sequence and adds the sequence and SMILES key into .csv file
triCenter = configBools[7]  #Centers in protein at trigram location instead of isolating trigrams for binding
NDMscreen = configBools[8] #Only extracts ligands from respective spreadsheet and dock against a given protein
    
#Variable Prep - initialies variables and path to Chembind excel sheet
mg = mygene.MyGeneInfo()  #instantiates mygene library
wb = xlrd.open_workbook(filename)
sheet = wb.sheet_by_index(0)
rows = sheet.nrows #Rows in sheet of protein + ligands
protnames = [None] * (rows  - 1)#array of protein HGNCs
ligSMILES = [None] * (rows  - 1) #array of ligand identifiers (SMILES)
lignames = [None] * (rows  - 1) #array of ligand names
ligCID = [None] * (rows  - 1) #array of ligand keys
colloqNames = [None] * (rows - 1) #Array of colloquial protein names

if CSVprepare == True: #Initializes needed settings if AIBind used to identify binding locations
    import csv
    CSVpath = '/home/' + user + '/Proteins/AIBindCSV.csv'
    ligChiKey = [None] * (rows - 1)
    
pdbqtNames = [] #For .pdbqt mmCIF subsections
pastfilename = '' #sets up pastfilename variable to stop unecessary downloading of repeat files
Protpath = '/home/' + user + '/Proteins'  #Paths to all major file locations, user-dependent
Deskpath = '/home/' + user
ligpath = '/home/' + user + '/Ligands'
indexpath = '/home/' + user + '/index.txt'
vinapath = '/home/' + user + '/vin.txt'
typepath = '/home/' + user + '/type.txt'
mmcifpath = '/home/' + user + '/mmcif.txt'
outpath = '/home/' + user + '/OutPoses'
endindex = rows - 1  #Used as a variable for for loops based off of chem spreadsheet
pfilename = ''  #Global variable for working with the protein file
uniprot = '' #uniprot ID for niche cases
p = '' # pure pdb ID without .pdb
alphaflag = False #used for proper .pdbqt processing
url = 'https://www.uniprot.org/uploadlists/' #For niche uniprot ID conversions
mygenefailed = False #For protein downloading from various sources
downloadfailed = False #For aborting the docking programs for proteins that do not exist
alphafold = 0 #Counts # of alphafold downloads
PDBdown = 0 #Counts # of pdb downloads
faileddownloads = 0 #Counts number of failed downloads
mmCIFflag = False  #Flags if downloaded file is of the mmCIF filetype
pdbTooLarge = False #Determines if pdb file needs to be broken up
mapfileName = 'All_Targets_DrugBank_BindingDB_DTC.xls'

#Function definitions - Functions used in main code runtime for processing and docking substrates
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
        
def fillarrays():
    global SMILE, ligID
    if NDMscreen == False:
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
                smiIndex = smiIndex + 1

            CIDFound = False
            CIDindex = 0
            try:
                while CIDFound == False:
                    if 'Chi' in sheet.cell_value(0, CIDindex) or 'Stitch' in sheet.cell_value(0, CIDindex):
                        ligCID[i] = sheet.cell_value(i + 1, CIDindex)
                        CIDFound = True
                        ligID = ligCID[i]
                    CIDindex = CIDindex + 1
            except:
                ligCID[i] = str(i)
                ligID = ligCID[i]
            
            if CSVprepare == True: #Gets original inChiKey for AIBind if so configured
                ligChiKey[i] = ligID
            ligCID[i] = (ligCID[i][4:]).lstrip('0') #Processes ligand ID for pcp
        else:
            ligSMILES[i] = sheet.cell_value(i + 1, 3)
            ligCID[i] = sheet.cell_value(i + 1, 1)
            colloqNames[i] = sheet.cell_value(i + 1, 4)
            lignames[i] = sheet.cell_value(i + 1, 2)
            ligID = ligCID[i]
    else:
        smiFound = False
        smiIndex = 0
        while smiFound == False:
            if 'smi' in sheet.cell_value(0, smiIndex):
                ligSMILES[i] = sheet.cell_value(i + 1, smiIndex)
                smiFound = True
            smiIndex = smiIndex + 1
        
        chiKeyFound = False
        chiKeyIndex = 0
        while chiKeyFound == False:
            if 'ick' in sheet.cell_value(0, chiKeyIndex):
                ligCID[i] = sheet.cell_value(i + 1, chiKeyIndex)
                chiKeyFound = True
                ligID = ligCID[i]
            chiKeyIndex = chiKeyIndex + 1
            
def downloadprotein(i): #Converts HGNC to uniprot name, downloads PDB file from uniprot
    global p, pfilename, uniprot, mmCIFflag, mygenefailed, alphaflag, downloadfailed, model1flag, alphafold, PDBdown, faileddownloads, pdbTooLarge, recpdbqt
    mygenefailed = False #resets indexer at start of function
    downloadfailed = False
    pdbFound = False
    model1flag = False #Lets the program know that the alphaprot file is named model 1
    os.chdir(Protpath)

    if absnegativeFormat == True:
        
        mapfilename = os.path.join(here, mapfileName)
        workbook = xlrd.open_workbook(mapfilename)
        sheets = workbook.sheet_by_index(0)
        srows = sheets.nrows
        for b in range(0, srows):
            if colloqNames[i] in str(sheets.cell_value(b + 1, 3)):
                uniprotID = sheets.cell_value(b + 1, 4)
                break
        out = mg.querymany(uniprotID, scopes = 'uniprot', fields = 'pdb', species = 'human')
    else:      
        HGNC = protnames[i]  #ID for genecard symbol from excel sheet
        out = mg.querymany(HGNC, scopes = 'symbol', fields = 'pdb, uniprot', species = 'human') #Checks the mygene archive
        #print(out) #for verification purposes
    
    try:
        protein = out[0]['pdb']   #Attempts to read PDB ID from uniprot results
        pdbFound = True
    except KeyError:
        mygenefailed = True
        if absnegativeFormat == False:
            try:
                uniprot = out[0]['uniprot']['Swiss-Prot']
            except KeyError:  #Uses uniprot library to convert Genecard ID to Uniprot Accession ID
                try:
                    os.system('scrapy fetch https://legacy.uniprot.org/uniprot/?query=gene:' + HGNC + ' > out.txt')
                            
                    with open('out.txt', 'r') as f:
                        Data = f.read() 

                    firstResultIndex = Data.index('selectAll-resultStatus')
                    uniprot = Data[firstResultIndex + 58:firstResultIndex + 64]

                    os.remove('out.txt')
                    with open('slub.txt', 'w') as f:
                        f.write('CUUM')
                    mygenefailed = True
                except:
                    downloadfailed = True
        else:
            uniprot = uniprotID
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
                        prot = out2[0]['uniprot']['Swiss-Prot']
                    else:
                        prot = uniprotID
                   
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    alphaname = 'AF-' + prot + '-F1-model_v2.pdb'
                    request.urlretrieve(url + alphaname, alphaname)
                    os.rename(alphaname, prot + '.pdb')
                    pfilename = prot + '.pdb'
                    p = prot
                    alphaflag = True
            else:
                #Downloads predicted alphafold structure if PDB ID not recognized
                if absnegativeFormat == False:
                    out2 = mg.querymany(HGNC, scopes = 'symbol', fields = 'uniprot', species = 'human')
                    prot = out2[0]['uniprot']['Swiss-Prot']
                else:
                    prot = uniprotID
               
                url = 'https://alphafold.ebi.ac.uk/files/'
                alphaname = 'AF-' + prot + '-F1-model_v2.pdb'
                request.urlretrieve(url + alphaname, alphaname)
                os.rename(alphaname, prot + '.pdb')
                pfilename = prot + '.pdb'
                p = prot
                alphaflag = True
    else:  #Downloads predicted alphafold structure if no PDB ID found
        try:
            if absnegativeFormat == False:
                pfilename = uniprot + '.pdb'
            else:
                uniprot = uniprotID
               
            url = 'https://alphafold.ebi.ac.uk/files/'
            alphaname = 'AF-' + uniprot + '-F1-model_v2.pdb'
            request.urlretrieve(url + alphaname, alphaname)
            os.rename(alphaname, uniprot + '.pdb')
            pfilename = uniprot + '.pdb'
            p = uniprot
            alphaflag = True
        except:
            print("Protein " + HGNC + " not found in any database!")
            downloadfailed = True
   
       
class breakLoop(Exception): pass

def receptorprep(pdb, pdbid, index):#Prepares protein receptor for docking, can set up code reading downloaded PDB name for pdbid
    global pfilename, subsection2Array, pdbqtNames, pdbname, recpdbqt, diffx, diffy, diffz, firstWrite, pdbTooLarge, residueNumber
    pdbname = pdb
    indexnum = 0
    overflow = 0
    diffx = 0
    diffy = 0
    diffz = 0
    os.chdir(Protpath) #Converts receptor to PDBQT file
    recpdbqt = str(pdbid + '.pdbqt')

    with open(pdb, 'r') as f:
            
        size = len(f.readlines())
        try:
            if data != 'SUCCESS':
                if size > 100000:  #File size where PDB protein is split up into subsections
                    pdbTooLarge = True
                elif testing != 'SUCCESS':
                    pdbTooLarge = True
                else:
                    pdbTooLarge = False
        except:
            True
            #Residue code here
           
    if mmCIFflag == True and residueBreak == False:   #Breaks up mmCIF file into subsections
        mmCIFtype = pdbid + '.pdb'
        myFile = open(mmCIFtype, 'r')
        lines = myFile.readlines()
        linelen = len(lines)
       
        subsectionArray = [] #Array of mmCIF subsections
        subsectionLength = 50000 #Multiplier for length of each subsection (3  = 30,000 lines)
        overlap = 1000 #Multiplier for overlap of subsections (2 = 2000 lines {~1000 atoms})
    
        sectionNum = int(linelen / subsectionLength) + 1
        for secIndex in range(0, sectionNum):
            subsectionArray.append([])
        
        firstWrite = False
        d = 0
            
        def reformat(num):   #Main function for rewriting PDB numbers
            global indexer, firstWrite
            if firstWrite == False:
                indexer = 0
                firstWrite = True
            if d < subsectionLength * (num + 1):
                if lines[d][0:4] == 'ATOM' or lines[d][0:6] == 'HETATM':
                    indexer = indexer + 1
                    if indexer < 100000:  #Requires rewriting file format for proper conversion after atom 100,000
                        space = '     '
                        space = space[:-len(str(indexer))]
                        lines[d] = lines[d].replace(lines[d][6:12], space + str(indexer) + ' ')
                        subsectionArray[num].append(lines[d])
                    else:
                        space = '     '
                        space = space[:-len(str(indexer))]
                        lines[d] = lines[d].replace(lines[d][6:12], space + str(indexer))
                        subsectionArray[num].append(lines[d])
                else:
                    subsectionArray[num].append(lines[d])
            else:
                raise breakLoop
   
               
        subsecNum = 0
        while d < linelen:
            try:
                while d < linelen:
                    reformat(subsecNum)
                    d = d + 1
            except breakLoop:
                subsecNum = subsecNum + 1
                d = d - overlap
       
               
        subsection2Array = [x for x in subsectionArray if x != []]
        subsection2Array.sort(key = len)
    if mmCIFflag == False or residueBreak == True: 
        if pdbTooLarge == True or residueBreak == True:   #Processes the file to pdbqt format ahead of subsectioning
            residueNumber = 0
            if data != 'FAILED':  #In the future, modify to testing
                
                recpdbqt = pdbid + '.pdbqt'

                if mmCIFflag == True:
                    mmCIFtype = pdbid + '.pdb'
                    myFile = open(mmCIFtype, 'r')
                    lines = myFile.readlines()
                    linelen = len(lines)
                   
                    subsectionArray = [] #Array of mmCIF subsections
                    subsectionLength = 50000 #Multiplier for length of each subsection (3  = 30,000 lines)
                    overlap = 1000 #Multiplier for overlap of subsections (2 = 2000 lines {~1000 atoms})
                
                    sectionNum = int(linelen / subsectionLength) + 1
                    for secIndex in range(0, sectionNum):
                        subsectionArray.append([])
                    
                    firstWrite = False
                    d = 0
                        
                    def reformat(num):   #Main function for rewriting PDB numbers
                        global indexer, firstWrite
                        if firstWrite == False:
                            indexer = 0
                            firstWrite = True
                        if d < subsectionLength * (num + 1):
                            if lines[d][0:4] == 'ATOM' or lines[d][0:6] == 'HETATM':
                                indexer = indexer + 1
                                if indexer < 100000:  #Requires rewriting file format for proper conversion after atom 100,000
                                    space = '     '
                                    space = space[:-len(str(indexer))]
                                    lines[d] = lines[d].replace(lines[d][6:12], space + str(indexer) + ' ')
                                    subsectionArray[num].append(lines[d])
                                else:
                                    space = '     '
                                    space = space[:-len(str(indexer))]
                                    lines[d] = lines[d].replace(lines[d][6:12], space + str(indexer))
                                    subsectionArray[num].append(lines[d])
                            #else:
                                #subsectionArray[num].append(lines[d])
                        else:
                            raise breakLoop
               
                           
                    subsecNum = 0
                    while d < linelen:
                        try:
                            while d < linelen:
                                reformat(subsecNum)
                                d = d + 1
                        except breakLoop:
                            subsecNum = subsecNum + 1
                            d = d - overlap
                   
                           
                    subsection2Array = [x for x in subsectionArray if x != []]
                    with open(pdbid + '.pdbqt', 'w') as a:
                        for item in subsection2Array:
                            for line in item:
                                a.write(line)
                else:
                    try:    
                        pdbqtconv = 'prepare_receptor4.py -r ' + pfilename + ' ' + recpdbqt + ' -A hydrogens'
                        os.system(pdbqtconv)
                    except:
                        os.system('obabel -.pdb ' + pfilename + ' -opdbqt ' + recpdbqt)
                try:
                    os.rename(pdbid + '_model1.pdbqt', pdbid + '.pdbqt')
                except:
                    True
                with open(recpdbqt, 'r') as hydroTest:
                    hydroData = hydroTest.readlines()
                    

                hydroIndexes = []                
                with open(recpdbqt, 'w') as hydroWrite:
                    num = 0
                    for line in hydroData:
                        if line[0:4] == 'ATOM' and ' ' not in line[13:21]:  #Removes space from polar hydrogen coordindates for correct processing
                            clist = list(map(lambda i:i, line))
                            del clist[29]
                            string = ''.join([str(item) for item in clist])
                            hydroWrite.write(string)
                            hydroIndexes.append(num)
                        elif HETremove == True and line[0:6] == 'HETATM':  #Enables removal of HETATMs if code configured to
                            del line
                        else:
                            hydroWrite.write(line)
                        num = num + 1
                    
                residueSize = []
                if residueBreak == True:
                    
                    with open(recpdbqt, 'r') as f:
                        residueData = f.readlines()
                        residueAtoms = [[]]

                    count = 0 
                    resCount = 0
                    resName = ''
                    pastResName = ''
                    for resIndex in range(0, len(residueData)):   #Finds residue name from ATOM string
                        if resCount not in hydroIndexes:
                            resName = residueData[resIndex][17:21]
                        else:
                            resName = residueData[resIndex][18:22]
                                
                        if count == 0:
                            pastResName = resName

                        if pastResName == resName:
                            residueAtoms[resCount].append(residueData[resIndex])
                        else:
                            residueSize.append(len(residueAtoms[resCount]))
                            resCount = resCount + 1
                            residueAtoms.append([])
                            residueAtoms[resCount].append(residueData[resIndex])
                        
                        pastResName = resName
                        count = count + 1
                    
                    residueNumber = len(residueSize)
                    
                    if trigramCombine == True:  #Combines residues together into trigrams
                        trigrams = []
                        triSize = []
                            
                        for z in range(0, int(residueNumber / 3)):
                            trigrams.append([])
                            for r in range(0, 3):
                                for line in residueAtoms[(z * 3) + r]:
                                    trigrams[z].append(line)
                            triSize.append(len(trigrams[z]))
                            
                        trigramNumber = len(triSize)
                    
            else:
                recpdbqt = pdb

            
            if residueBreak == False:
                if testing != 'FAIL':
                    subsectionLength = 20000 #Multiplier for length of each subsection
                    overlap = 1000
                else:
                    subsectionLength = 10000
                    overlap = 500
                    
                with open(recpdbqt, 'r') as f:
                    lines = f.readlines()
                    linelen = len(lines)
                        
                subsectionArray = []
                sectionNum = int(linelen / subsectionLength) + 1 #Automatically generates array of subsections
                for secIndex in range(0, sectionNum):
                    subsectionArray.append([])

            else:  #Correctly prepares length of subsectionArray in case of residue/trigram breakage
                subsectionArray = []
                if trigramCombine == False:
                    for secIndex in range(0, residueNumber):
                        subsectionArray.append([])
                else:
                    for secIndex in range(0, trigramNumber):
                        subsectionArray.append([])  

            
            firstWrite = False
            d = 0

            subsecNum = 0
            def reformat(num):   #Main function for rewriting PDB numbers to subsections
                global indexer, firstWrite
                if firstWrite == False:
                    indexer = 0
                    firstWrite = True
                if residueBreak == False:
                    if d < subsectionLength * (num + 1):
                        if lines[d][0:4] == 'ATOM' or lines[d][0:6] == 'HETATM':
                            indexer = indexer + 1
                            space = '     '
                            space = space[:-len(str(indexer))]
                            lines[d] = lines[d].replace(lines[d][6:12], space + str(indexer) + ' ')
                            subsectionArray[num].append(lines[d])
                        else:
                            subsectionArray[num].append(lines[d])
                    else:
                        raise breakLoop
                else:
                    if trigramCombine == False:
                        if d < residueSize[num]:
                            if residueAtoms[num][d][0:4] == 'ATOM' or residueAtoms[num][d][0:6] == 'HETATM':
                                indexer = indexer + 1
                                space = '     '
                                space = space[:-len(str(indexer))]
                                atom = residueAtoms[num][d].replace(residueAtoms[num][d][6:12], space + str(indexer) + ' ')
                                subsectionArray[num].append(atom)
                            else:
                                subsectionArray[num].append(residueAtoms[num][d])
                        else:
                            raise breakLoop
                    else:
                        if d < triSize[num]:
                            if trigrams[num][d][0:4] == 'ATOM' or trigrams[num][d][0:6] == 'HETATM':
                                indexer = indexer + 1
                                space = '     '
                                space = space[:-len(str(indexer))]
                                atom = trigrams[num][d].replace(trigrams[num][d][6:12], space + str(indexer) + ' ')
                                subsectionArray[num].append(atom)
                            else:
                                subsectionArray[num].append(trigrams[num][d])
                        else:
                            raise breakLoop
                
            if residueBreak == False:
                while d < linelen:
                    try:
                        while d < linelen:
                            reformat(subsecNum)
                            d = d + 1
                    except breakLoop:
                        subsecNum = subsecNum + 1
                        d = d - overlap
                
            else:
                if trigramCombine == False:
                    while subsecNum < residueNumber:
                        try:
                            while True:
                                reformat(subsecNum)
                                d = d + 1
                        except breakLoop:
                            d = 0
                            subsecNum = subsecNum + 1
                else:
                    while subsecNum < trigramNumber:
                        try:
                            while True:
                                reformat(subsecNum)
                                d = d + 1
                        except breakLoop:
                            d = 0
                            subsecNum = subsecNum + 1
                            
            subsection2Array = [x for x in subsectionArray if x != []]  
            
            if data == 'FAILED':
                d = 0
                pdbqtNames = []
                pdbNames = []
                for item in subsection2Array:
                    d = d + 1
                    pdbqtName = str(d) + '.pdbqt'
                    fName = str(d) + '.pdb'
                    pdbqtNames.append(pdbqtName)
                    pdbNames.append(fName)
                    with open(fName, 'w') as file:
                        for line in item:
                            file.write(line)
                        
                        pdbqtconv = 'prepare_receptor4.py -r ' + fName + ' ' + pdbqtName + ' -A hydrogens'

                    try:
                        os.system(pdbqtconv)
                    except:
                        os.system('obabel -.pdb ' + pfilename + ' -opdbqt ' + recpdbqt)
                    try:
                        os.rename(str(d) + '_model1.pdbqt', str(d) + '.pdbqt')
                    except:
                        pass

            else:
                d = 0   #Writes .pdbqt files for each subsection
                pdbqtNames = []
                for item in subsection2Array:
                    d = d + 1
                    pdbqtName = str(d) + '.pdbqt'
                    pdbqtNames.append(pdbqtName)
                    with open(pdbqtName, 'w') as file:
                        for line in item:
                            file.write(line)
               
            #os.remove(recpdbqt)
        else:
            recpdbqt = pdbid + '.pdbqt'
            try:
                pdbqtconv = 'prepare_receptor4.py -r ' + pdb + ' ' + recpdbqt + ' -A hydrogens'
                os.system(pdbqtconv)
                
            except:
                os.system('obabel -.pdb ' + pfilename + ' -opdbqt ' + recpdbqt)
            try:
                os.rename(pdbid + '_model1.pdbqt', recpdbqt)
            except:
                pass
                
            with open(recpdbqt, 'r') as hydroTest:
                hydroData = hydroTest.readlines()
            with open(recpdbqt, 'w') as hydroWrite:
                for line in hydroData:
                    if line[0:4] == 'ATOM' and ' ' not in line[13:21]:
                        clist = list(map(lambda i:i, line))
                        del clist[29]
                        string = ''.join([str(item) for item in clist])
                        hydroWrite.write(string)
                    else:
                        hydroWrite.write(line)
                
            if data == 'FAILED':
                newData = []
                with open(recpdbqt, 'r') as f:
                    dataList = f.readlines()
                    for line in dataList:
                        if line[0:6] == 'HETATM':
                            atomType = line[11:15].strip(' ')
                            if atomType == 'AS':
                                continue
                        newData.append(line)
                with open(recpdbqt, 'w') as f:
                    for line in newData:
                        if HETremove == True:
                            if line[0:6] != 'HETATM':
                                f.write(line)
                        else:
                            f.write(line)
   
                               
    else:
        d = 0   #pdbqt conversion for mmCIF subsections
        pdbNames = []
        for item in subsection2Array:
            d = d + 1
            fName = str(d) + '.pdb'
            pdbqtName = str(d) + '.pdbqt'
            pdbqtNames.append(pdbqtName)
            pdbNames.append(fName)
            with open(fName, 'w') as f:
                for line in item:
                    f.write(line)
                if d < len(subsection2Array):
                    f.write('END')
            pdbqtconv = 'prepare_receptor4.py -r ' + fName + ' ' + pdbqtName + ' -A hydrogens'
            try:
                os.system(pdbqtconv)
            except:
                os.system('obabel -.pdb ' + pfilename + ' -opdbqt ' + recpdbqt)
           
    #if residueBreak == False and pdbTooLarge == False and mmCIFflag == False:
    #    with open(recpdbqt, 'r') as r:
    #        fileData = r.readlines()
    #    with open(recpdbqt, 'w') as r:
    #        for line in fileData:
    #           if HETremove == False:
    #                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
    #                    r.write(line)
    #            else:
    #                if line[0:4] == 'ATOM':
    #                  r.write(line)
                    

#################################### AIBIND INTEGRATION SECTION ################################################
    if CSVprepare == True and residueBreak == True: #Identifies binding locations on protein
        os.chdir(Protpath)
        with open('bindingSites.txt', 'r') as f:
            AIBResults = f.read()
            
        AIBResults = ast.literal_eval(AIBResults)
        with open('dub.txt', 'w') as f:
            f.write(str(index - docking_range[0]) + '\n')
            f.write(str(AIBResults))
        AIBResults = AIBResults[index - docking_range[0]]
        #with open('dub.txt', 'w') as f:
         #   f.write(str(AIBResults))
            
        #try:
        #    AIBResults = AIBResults[index - docking_range[0]]
        #except:
        #    with open('ERROR.txt') as f:
        #        print('INDEXING ERROR!')
        

        trigramPockets = []
        indexer = 0
        Overlap = 1
        for num in range(0, len(AIBResults)): #for each set of pockets
            if len(AIBResults[num]) > 0:
                trigramPockets.append([])
                for c in range(AIBResults[num][0] - Overlap, AIBResults[num][1] + 1 + Overlap):
                    if c < len(subsectionArray):
                        trigramPockets[num].append(subsectionArray[c])
            
        pdbqtIndexer = 0
        for item in trigramPockets:
            if len(item) == 0:
                del trigramPockets[pdbqtIndexer]
            pdbqtIndexer = pdbqtIndexer + 1
    
        pdbqtNames = []
        if triCenter == False:
            global triSections
            for c in range(0, len(trigramPockets)): #rewrites binding pockets as pdbqt files
                with open('trigram' + str(c) + '.pdbqt', 'w') as f:
                    for item in trigramPockets[c]:
                        for line in item:
                            f.write(line)
                pdbqtNames.append('trigram' + str(c) + '.pdbqt')
        elif triCenter == True:
            triSections = []
            for c in range(0, len(trigramPockets)): #rewrites binding pockets as pdbqt files
                with open('trigram' + str(c) + '.pdbqt', 'w') as f:
                    for item in trigramPockets[c]:
                        for line in item:
                            f.write(line)
                triSections.append('trigram' + str(c) + '.pdbqt')    
                pdbqtNames.append(recpdbqt)
    
        if len(pdbqtNames) == 0:
            with open('noPockets.txt', 'w') as f:
                f.write('temp')
    
def downloadligand(): #Downloads ligand file
    if SMILEUse == False and NDMscreen == False:
        os.chdir(ligpath)
        try:
            pcp.download('SDF', ligandfilename , int(ligCID[i]), 'cid', overwrite = True)
        except:
            print("Ligand not found in PubChem database!")
    lipath = ligpath + '/'
    lipath = lipath[:-1]


def ligandprep(): #Prepares ligand for docking
    global lipdbqtname
    os.chdir(ligpath)

    lipdbqtname = str(ligCID[i] + '.pdbqt')
    
    removepath = 'rm '
    if SMILEUse == False and NMDScreen == False:
        lipdbname = 'temp'
        os.system('obabel -isdf ' + ligandfilename + ' -O ' + ligandfilename + ' -m --gen3d best -p') #NOTE - for ADT may need to add Gasteiger
        os.system('obabel ' + ligandilename + ' -O ' + lipdbname + '.pdb -m')
        os.system('prepare_ligand4.py -l temp1.pdb -o ' + lipdbqtname) #-A hydrogens bonds
    else:
        try:
            #lipdbname = 'temp'
            pcp.download('SDF', ligandfilename, ligCID[i], 'inchikey', overwrite = True)
            
            os.system('obabel ' + ligandfilename + ' -O ' + lipdbqtname + ' -m') 
            #os.system('prepare_ligand4.py -l ' + lipdbname + ' ' + lipdbqtname + ' -A hydrogens')
            try:
                f = open(lipdbqtname, 'r')
            except:
                lipdbqtname = ligCID[i] + '1.pdbqt'
                #os.rename('temp.pdbqt', lipdbqtname)
            
        except:
            #if lipdbqtname == '.pdbqt':
    
            with open('SMILE.txt', 'w') as f:
                f.write(ligSMILES[i])
            os.system('obabel -ismi SMILE.txt -O SMILE.pdbqt -m --gen3d best -p')  
            #os.system('prepare_ligand4.py -l ' + lipdbname + ' ' + lipdbqtname + ' -A hydrogens') 
            #os.system('obabel -ismi SMILE.txt -opdbqt -m --gen3d best -p -add hydrogens')
            try:
                os.rename('SMILE.pdbqt', lipdbqtname)
            except:
                try:
                    os.rename('SMILE1.pdbqt', lipdbqtname)
                except:
                
                    pass
                  
            removepath = removepath + 'SMILE.pdbqt'
            #os.system(removepath)
           
    #sys.exit()
    
class coordExtract: #Master superclass for extracting coordinate data
    xcoords = []
    ycoords = []
    zcoords = []
    coordData = ''
    lineNum = 0
    xcenter = 0
    ycenter = 0
    zcenter = 0
    xsize = 0
    ysize = 0
    zsize = 0
    def __init__(self, fileName, fileType):
        self.fileName = fileName
        self.fileType = fileType
        
    def readData(self):
        tempData = []
        with open(self.fileName, 'r') as f:
            tempData = f.readlines()
        for c in range(0, len(tempData)):
            if tempData[c][0:4] != 'ATOM' and tempData[c][0:6] != 'HETATM':
                tempData[c] = ''
        self.coordData = ''.join(tempData)
    
    def coordExtract(self):
        self.lineNum = int((len(self.coordData) - 28) / 80)
        for i in range(0, self.lineNum):
            self.xcoords.append(float(self.coordData[30 + i * 80:38 + i * 80]))
            self.ycoords.append(float(self.coordData[39 + i * 80:46 + i * 80]))
            self.zcoords.append(float(self.coordData[47 + i * 80:54 + i * 80]))

        
    def calcCenter(self):
        self.xcenter = statistics.mean(self.xcoords)
        self.ycenter = statistics.mean(self.ycoords)
        self.zcenter = statistics.mean(self.zcoords)

    def calcGridSize(self):
        self.xsize = max(self.xcoords) - min(self.xcoords) 
        self.ysize = max(self.ycoords) - min(self.ycoords) 
        self.zsize = max(self.zcoords) - min(self.zcoords) 

class ligCoordExtract(coordExtract):
    pass

class Dock:
    def __init__(self, lipdbqtname, Protname, xcenter, ycenter, zcenter, xsize, ysize, zsize, optimize): #optimize as bool
        self.lipdbqtname = lipdbqtname
        self.Protname = Protname
        self.xcenter = xcenter
        self.ycenter = ycenter
        self.zcenter = zcenter
        self.xsize = xsize
        self.ysize = ysize
        self.zsize = zsize
        self.optimize = optimize

    def Dock(self):
        v = Vina(sf_name = 'vina')
        os.chdir(Protpath)
        v.set_receptor(self.Protname)
        os.chdir(ligpath)
        v.set_ligand_from_file(self.lipdbqtname)
        v.compute_vina_maps(center = [self.xcenter, self.ycenter, self.zcenter], box_size = [self.xsize, self.ysize, self.zsize])
        os.chdir(Deskpath)
        if self.optimize == True: #Test with ligand optimization
            if data != 'FAILED':
                lowest = v.optimize()
                v.dock(exhaustiveness = 20, n_poses = 20)
            else:
                if testing == 'SUCCESS':
                    v.dock(exhaustiveness = 16, n_poses = 20)
                else:
                    v.dock(exhaustiveness = 8, n_poses = 20)
        else:
            v.dock(exhaustiveness = 16, n_poses = 20)
            
        #os.chdir(outpath)
        #if CSVprepare == False:
            #v.write_poses(outposeName, n_poses = 1, overwrite = True)
        #else:
            #v.write_poses(outposeName[dockCount], n_poses = 1, overwrite = True)
                                
def autoDock(protein, ligand):
    os.chdir(ligpath)
    Ligand = ligCoordExtract(ligand, 'ligand')
    Ligand.readData()
    Ligand.coordExtract()
    Ligand.calcGridSize()
    
    #with open('ligSizeCheck.txt', 'w') as f: #For testing
        #f.write(str(Ligand.xsize))
        
    os.chdir(Protpath)
    if isinstance(protein, list) == False:
        Protein = coordExtract(protein, 'protein')
        Protein.readData()

        Protein.xcoords = []
        Protein.ycoords = []
        Protein.zcoords = []
            
        Protein.coordExtract()
        Protein.calcCenter()
        Protein.calcGridSize()
    
        
        xsize = Protein.xsize + Ligand.xsize + 6
        ysize = Protein.ysize + Ligand.ysize + 6
        zsize = Protein.zsize + Ligand.zsize + 6
        xcenter = Protein.xcenter
        ycenter = Protein.ycenter
        zcenter = Protein.zcenter

        with open('coordCheck.txt', 'w') as c: #For testing
            c.write(str(xsize) + ' ' + str(ysize) + ' ' + str(zsize))
            c.write(str(xcenter) + ' ' + str(ycenter) + ' ' + str(zcenter))
            
        soloDock = Dock(Ligand.fileName, Protein.fileName, \
                    Protein.xcenter, Protein.ycenter, Protein.zcenter, \
                        xsize, ysize, zsize, False)
        soloDock.Dock()

    else:
        protLen = len(protein)
        protclassList = []
        xcenter = []
        ycenter = []
        zcenter = []
        xsize = []
        ysize = []
        zsize = []

        global dockCount
        dockCount = 0
        with open('trigramNum.txt', 'w') as f:
            f.write(str(protLen))
        for item in protein:
            os.chdir(Protpath)
            Protein = coordExtract(item, 'protein')
            
            Protein.readData()

            Protein.xcoords = []
            Protein.ycoords = []
            Protein.zcoords = []

            Protein.xcoords = []
            Protein.coordExtract()

            try:  
                Protein.calcCenter()
            except:
                break
            Protein.calcGridSize()

            with open('coordCheck.txt', 'w') as c: #For testing
                c.write(str(Protein.xsize) + '\n')
                c.write(str(Protein.ysize) + '\n')
                c.write(str(Protein.zsize) + '\n')

            if triCenter == False:
                xsize.append(Protein.xsize + Ligand.xsize + 6)
                ysize.append(Protein.ysize + Ligand.ysize + 6)
                zsize.append(Protein.zsize + Ligand.zsize + 6)
            else:
                xsize.append(Protein.xsize + 3)
                ysize.append(Protein.ysize + 3)
                zsize.append(Protein.zsize + 3)
            xcenter.append(Protein.xcenter)
            ycenter.append(Protein.ycenter)
            zcenter.append(Protein.zcenter)
            

            #with open('Check.txt', 'a') as c: #For testing
                #c.write(Protein.fileName + '\n')
                #c.write(str(len(subsection2Array[0])) + '\n')
                
                
                #c.write(str(xsize[dockCount]) + ' ' + str(ysize[dockCount]) + ' ' + str(zsize[dockCount]) + '\n')
                #c.write(str(xcenter[dockCount]) + ' ' + str(ycenter[dockCount]) + ' ' + str(zcenter[dockCount]) + '\n')

            
            #multiDock = Dock(Ligand.fileName, Protein.fileName, \
                    #Protein.xcenter, Protein.ycenter, Protein.zcenter, \
                        #Protein.xsize, Protein.ysize, Protein.zsize, False)
            multiDock = Dock(Ligand.fileName, Protein.fileName, \
                    xcenter[dockCount], ycenter[dockCount], zcenter[dockCount], \
                        xsize[dockCount], ysize[dockCount], zsize[dockCount], False)
            multiDock.Dock()

            dockCount = dockCount + 1
            
        with open(mmcifpath, 'w') as file:
            file.write(str(len(pdbqtNames)))
            
            #os.chdir(Protpath)          #Writes output pose as a .pqbqt file
    #v.write_poses(pdbqt_filename="docking_results.pdbqt")
    #v.write_poses(lipdbqtname, n_poses = 1, overwrite = True)
    
def outfileCalc():
    global outposeName
    if CSVprepare == False: #given protein, modify for multiple proteins
        outposeName = ligCID[i] + '%' + p + '.pdbqt'
    else:
        outposeName = []
        indexer = 0
        for item in pdbqtNames:
            outposeName.append(ligID + '%' + pdbqtNames[indexer] + str(i) + '.pdbqt')
            indexer = indexer + 1
    
###############################  MAIN RUNTIME LOOP ############################################################################


        
with open(indexpath, 'r') as myFile:  #Reads current index and runs functions with corrosponding ligand + PDB pairs
    with open(vinapath, 'r') as f:   #Asseses whether program is being rerun due to failure
        data = f.read()
        
    try:   #Asseses whether program should be run in safe mode
        with open('safeMode.txt', 'r') as f:
            testing = f.read()
    except:
        testing = 'SUCCESS'
    
    index = myFile.readlines()
    i = int(index[0])

    
    protConverted = False #For single-protein NDM purposes
    fillarrays()

    ligandfilename = ligCID[i] + '.sdf'
    if ligandfilename != pastfilename:
        downloadligand()
        ligandprep()
    if NDMscreen == False:
        if data != 'FAILED':
            
            downloadprotein(i)
        else:
            
            try:
                os.chdir(Protpath)
                with open('protname.txt', 'r') as f:
                    pfilename = str(f.readlines())
                    p = pfilename.rstrip(pfilename[-4])
                    with open(pfilename, 'r') as f:
                        pass
            except:
                downloadprotein(i)
                
        if mmCIFflag == True:
            if pfilename[-4] == '.pdb':
                mmCIFflag = False
            else:
                pass
            
        with open('protname.txt', 'w') as f:
            f.write(pfilename)

        
    else: #Executes code for obtaining protein target
        os.chdir(Protpath)
        try:
            with open('protname.txt', 'r') as f:
                pfilename = f.readlines()
                p = pfilename.rstrip(pfilename[-4])
            with open(pfilename, 'r') as f:
                pass
            protConverted = True
            
        except:
        
            pdbFound = False
            mmCIFflag = False
            alphaflag = False
            model1flag = False
            downloadfailed = False
            PDBdown  = 0
            with open('newprot.txt', 'r') as c:
                protein = c.read()
            
            
            out = mg.querymany(protein, scopes = 'uniprot', fields = 'pdb', species = 'human')
            try:
                protein = out[0]['pdb']   #Attempts to read PDB ID from uniprot results
                pdbFound = True
            except KeyError:
                try:
                    job_id = submit_id_mapping(
                        fromDB="UniProtKB_AC-ID", toDB="PDB", ids=[protein])
                    protein = get_id_mapping_results(job_id)
                    pdbFound = True
                except:
                    pass
            if pdbFound == True:
                if isinstance(protein, list):     #Processes PDB IDs, attempts to download from RCSB
                    regex = re.compile('[^a-zA-Z0-9]')
                    RCSBID = [None] * len(protein)
                    for j in range(0, len(RCSBID)):
                        protein[j] = regex.sub('', protein[j])
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
                except Exception as e:
                    print(e)
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
                            prot = protein
                           
                            url = 'https://alphafold.ebi.ac.uk/files/'
                            alphaname = 'AF-' + prot + '-F1-model_v2.pdb'
                            request.urlretrieve(url + alphaname, alphaname)
                            os.rename(alphaname, prot + '.pdb')
                            pfilename = prot + '.pdb'
                            p = prot
                            alphaflag = True
                    else:
                        #Downloads predicted alphafold structure if PDB ID not recognized
                        prot = protein
                        url = 'https://alphafold.ebi.ac.uk/files/'
                        alphaname = 'AF-' + prot + '-F1-model_v2.pdb'
                        request.urlretrieve(url + alphaname, alphaname)
                        os.rename(alphaname, prot + '.pdb')
                        pfilename = prot + '.pdb'
                        p = prot
                        alphaflag = True
            else:  #Downloads predicted alphafold structure if no PDB ID found
                
                try:
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    alphaname = 'AF-' + protein + '-F1-model_v2.pdb'
                        
                    request.urlretrieve(url + alphaname, alphaname)
                    os.rename(alphaname, protein + '.pdb')
                    pfilename = protein + '.pdb'
                    p = protein
                    alphaflag = True

                except:
                    print("Protein " + protein + " not found in any database!")
                    downloadfailed = True
       
            with open('protname.txt', 'w') as f:
                f.write(pfilename)

    
    
    if downloadfailed != True:
        if protConverted == False or CSVprepare == True:
           
            receptorprep(pfilename, p, i)
        else:
            
            recpdbqt = p + '.pdbqt'
        if mmCIFflag == False and pdbTooLarge == False:   #Tests if code is able to translate to .pdbqt through obabel, if not conversion is done automatically
            try:
                with open(recpdbqt, 'r') as myFile:
                    True
                pfilename = recpdbqt
            except:
                if data != 'SUCCESS':
                    receptorprep(pfilename, p, i)
            
        if data != 'SUCCESS':
            for item in pdbqtNames:
                with open(item, 'r') as hydroTest:
                    hydroData = hydroTest.readlines()
                        
                hydroIndexes = []                
                with open(item, 'w') as hydroWrite:
                    num = 0
                    for line in hydroData:
                        if line[0:4] == 'ATOM' and ' ' not in line[13:21]:  #Processes polar hydrogens correctly
                            clist = list(map(lambda i:i, line))
                            del clist[29]
                            string = ''.join([str(item) for item in clist])
                            hydroWrite.write(string)
                            hydroIndexes.append(num)
                        else:
                            hydroWrite.write(line)
                        num = num + 1

        outfileCalc()
        
        if pdbTooLarge == False and mmCIFflag == False and residueBreak == False:
           
            autoDock(pfilename, lipdbqtname)
        else:
            #with open('mmur.txt', 'w') as e:
               #e.write(str(pdbqtNames))
            autoDock(pdbqtNames, lipdbqtname)
            
        
        #sys.exit()   #For parameter testing purposes
                
        pastfilename = ligandfilename #Variables used to prevent multiple downloads of the same file
        
    else:
        print('Protein ' + protnames[i] + ' has not been found!')  #If no valid protein found, no docking can be run

if NDMscreen == False:
    with open(typepath, 'w') as f:
        if alphaflag == True:
            f.write('ALPHAFOLD')
        elif alphaflag == False and downloadfailed == False:
            f.write('PDB')
        else:
            f.write('FAILED')