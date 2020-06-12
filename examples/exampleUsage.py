from DecoID.DecoID import DecoID
libFile = "../databases/HMDB_experimental.db"
libFile = "none"
key = ""
useAuto = True
numCores = 10#int(sys.argv[2])
file = "../exampleData/Asp-Mal_1uM_5Da.mzML"#sys.argv[1]

usePeaks = True
DDA = True#False
massAcc = 10
fragThresh= 0.01
offset = 25
useIso = True
threshold = 0
lam = 10
if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    decID = DecoID(libFile, useAuto, numCores,api_key=key)
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions="../exampleData/peak_table.csv")
    decID.identifyUnknowns()
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold)
