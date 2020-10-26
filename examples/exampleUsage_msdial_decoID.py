from DecoID.DecoID import DecoID

#sets database to use
libFile = "../databases/HMDB_experimental.db"

#mzCloud key if necessary
key = "none"
mzCloudLib = "reference"

#number of parallel processes to use
numCores = 5

#filename of query MS/MS data
file = "../exampleData/IROA_P1-6_DIA_test_pos1.mzML"

#filename of peak list
peakFile = "../exampleData/IROA_p1-6_peak_table_pos_v3.csv"

#set parameters
usePeaks = True
DDA = False #data is DIA
massAcc = 10 #ppm tolerance
fragThresh= 0.01 #require non-zero dot product threshold
offset = .5 #half of isolation window width. Only for non-thermo data
useIso = True #use predicted M+1 isotopolgoue spectra
threshold = 0 #minimum dot product for reporting
lam = 50.0 #LASSO parameter
rtTol = float("inf") #retention time tolerance for database, inf means ignore RT
fragCutoff = 1000 #intensity threshold for MS/MS peaks


if __name__ == '__main__':
    #create DecoID object
    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key)

    #deconvolve with DecoID
    #read in data
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions=peakFile,frag_cutoff=fragCutoff)

    #search spectra
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)

    #deconvolve with MS-DIAL

    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key,label="_msdial")

    fn = "../exampleData/IROA_P1-6_DIA_test_pos1.txt"

    #read in MS-DIAL output
    decID.readMS_DIAL_data(fn,"Positive",massAcc,peakFile) #set polarity to Positive or Negative

    #directly search output against database using same peak list
    decID.searchSpectra("y",float("inf"),percentPeaks=fragThresh,iso=useIso,threshold=threshold,rtTol = rtTol)

    #combine the resulting taking the best scores from DecoID and MS-DIAL
    decID.combineResultsAutoAnnotate([fn.replace(".txt","_msdial"),file.replace(".mzML","")],fn.replace(".txt","_combined"),numHits = 3)