from DecoID.DecoID import DecoID
libFile = "../databases/mzCloud_reference.db"
key = open("../housekeeping/mzCloud_api_key.txt").readline().rstrip()
mzCloudLib = "reference"
numCores = 10#int(sys.argv[2])
file = "../exampleData/Asp-40uM_Mal1uM_1Da_NCE30_750ms_NCE50.mzML"#sys.argv[1]

usePeaks = True
DDA = True#False
massAcc = 25
fragThresh= 0.01
offset = 25
useIso = True
threshold = 0
lam = 0
rtTol = float("inf")
if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key)
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions="../exampleData/peak_table.csv",frag_cutoff=1000)
    decID.identifyUnknowns(rtTol=rtTol)
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)
