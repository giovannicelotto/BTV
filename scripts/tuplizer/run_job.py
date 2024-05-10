import pandas as pd
import glob, sys, re, os
import random
import subprocess
import time

def main(nFiles):
    # Define name of the process, folder for the files and xsections
    nanoPath="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06"
    flatPath="/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/btv_ntuples/TTToHadronic2024May06Tuple"
    nanoFileNames = glob.glob(nanoPath+"/**/*.root", recursive=True)
    flatFileNames = glob.glob(flatPath+"/**/*.root", recursive=True)
    
    if len(flatFileNames)==len(nanoFileNames):
        return
        
    time.sleep(3)
    
    nFiles = nFiles if nFiles != -1 else len(nanoFileNames)
    if nFiles > len(nanoFileNames) :
        nFiles = len(nanoFileNames)
    #nFiles to be done
    doneFiles = 0
    for nanoFileName in nanoFileNames:
        if doneFiles==nFiles:
            break
        fileNumber = re.search(r'\D(\d{1,4})\.\w+$', nanoFileName).group(1)
        print(fileNumber)
        filePattern = flatPath+"/**/"+"TTToH_"+fileNumber+".root"
        matching_files = glob.glob(filePattern, recursive=True)
        #print("Checking for ", flatPath+"/**/"+process+"_"+fileNumber+".root")

        if matching_files:
            print("TToH"+"_"+fileNumber+".root present. Skipped")
            continue
        print('nanoFile name', nanoFileName)
        print('fileNumber', fileNumber)
        print('flatPath', flatPath)
        subprocess.run(['sbatch', '-J', "TTToH"+"%d"%random.randint(1, 20), '/t3home/gcelotto/BTV/scripts/tuplizer/job.sh', nanoFileName, fileNumber, flatPath])
        doneFiles = doneFiles+1
    return 

if __name__ == "__main__":
    nFiles      = int(sys.argv[1]) if len(sys.argv) > 1 else -1
    main(nFiles)