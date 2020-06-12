from DecoID.DecoID import DecoID

libFile = "../databases/HMDB_experimental.db"

useAuto = True
numCores = 10#int(sys.argv[2])
file = "../exampleData/"#sys.argv[1]

usePeaks = True
DDA = True#False
massAcc = 10
fragThresh= 0.01
offset = 25
useIso = True
threshold = 0
lam = 10
if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")

    decID = DecoID.DecoID(False, libFile, useAuto, numCores, label=str(lam) + "_")
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions="../exampleData/NIST_1950/NIST1950_neg/peak_table.csv")
    decID.identifyUnknowns()
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold)
