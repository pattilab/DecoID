import sys
import os
#libFile = "none"
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src/DecoID/"))
from DecoID import DecoID
import multiprocessing



useAuto = False
numCores = int(sys.argv[3])#
file = sys.argv[1]#"../exampleData/Asp-40uM_Mal1uM_1Da_NCE30_750ms_NCE50.raw"
peakFile = sys.argv[4]
libFile = sys.argv[5]
usePeaks = True
DDA = False
massAcc = 10
useIso = False
numFiles = int(sys.argv[2])
if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    decID = DecoID.DecoID(False, libFile, useAuto, numCores)
    decID.readData(file, 2, usePeaks, DDA, massAcc,peakDefinitions=peakFile)
    #decID.identifyUnknowns(iso=useIso)
    DecoID.DecoID.prepareForCluster(decID, numFiles)

