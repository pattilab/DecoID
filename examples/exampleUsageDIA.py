from DecoID.DecoID import DecoID

#sets database to use
libFile = "../databases/HMDB_experimental.db"
libFile = "../databases/IROA_updated_NCE40.db"

#mzCloud key if necessary
key = "none"
mzCloudLib = "reference"

#number of parallel processes to use
numCores = 15

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
useIso = False #use predicted M+1 isotopolgoue spectra
threshold = 0 #minimum dot product for reporting
lam = 50.0 #LASSO parameter
rtTol = 0.5#float("inf") #retention time tolerance for database, inf means ignore RT
fragCutoff = 1000 #intensity threshold for MS/MS peaks

if __name__ == '__main__':

    #create DecoID object
    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key)

    #read in data
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions=peakFile,frag_cutoff=fragCutoff)

    #search spectra
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)
