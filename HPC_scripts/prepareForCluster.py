import sys
from DecoID import DecoID
import multiprocessing


useAuto = False
numCores = int(sys.argv[3])#
file = sys.argv[1]
peakFile = sys.argv[4]
libFile = sys.argv[5]
usePeaks = True
DDA = True
massAcc = 10
useIso = True
numFiles = int(sys.argv[2])
if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    decID = DecoID.DecoID(libFile, useAuto, numCores)
    decID.readData(file, 2, usePeaks, DDA, massAcc,peakDefinitions=peakFile)
    decID.identifyUnknowns(iso=useIso)
    DecoID.DecoID.prepareForCluster(decID, numFiles)

