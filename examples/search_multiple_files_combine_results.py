import DecoID.DecoID as DecoID
import numpy as np
import os

#parameters

#mzML file dir
datadir = "../exampleData/mzMLs/"

#db file path
libFile = "../databases/HMDB_experimental.db"

#number of parallel processes to use
numCores = 10

#filename of peak list
peakfile = "../exampleData/IROA_p1-6_peak_table_pos_v3.csv"

#set parameters
usePeaks = True #use MS1 data or not
DDA = True #data is DDA
massAcc = 10 #ppm tolerance
fragThresh= 0.01 #require non-zero dot product threshold
offset = .5 #half of isolation window width. Only for non-thermo data
useIso = True #use predicted M+1 isotopolgoue spectra
threshold = 0 #minimum dot product for reporting
lam = 5.0 #LASSO parameter
rtTol = float("inf") #retention time tolerance for database, inf means ignore RT
fragCutoff = 1000 #intensity threshold for MS/MS peaks
numHits = 3 #number of hits to keep for each feature

if __name__ == "__main__":
    
    decoID = DecoID.DecoID(libFile, "reference", numCores=numCores) #make DecoID object
    filenamesToCombine = [datadir + x.replace(".mzML","") for x in os.listdir(datadir) if ".mzML" in x] #get filenames to search with mzML ending removed
    for filename in [x for x in os.listdir(datadir) if ".mzML" in x]: #iterate through files
        print(datadir + filename) #print current filename
        decoID.readData(datadir + filename,2,usePeaks,DDA,massAcc,peakDefinitions=peakfile,frag_cutoff=fragCutoff) #read in data
        decoID.identifyUnknowns(resPenalty=lam,iso=useIso,rtTol=rtTol) #build on-the-fly unknown library
        decoID.searchSpectra("y",resPenalty=lam,iso=useIso,rtTol=rtTol) #deconvolve and search
    decoID.combineResultsAutoAnnotate(filenamesToCombine,datadir + "combined_resultsTOP" + str(numHits) + "_above" + str(threshold) + ".csv",numHits=numHits,min_score=threshold)