import sys
import os

#libFile = sys.argv[1]

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src/DecoID/"))
from DecoID import DecoID
import multiprocessing

filename = sys.argv[2]

useAuto = False
numCores = int(sys.argv[3])
peakFile = sys.argv[4]

usePeaks = bool(int(sys.argv[5]))
DDA = bool(int(sys.argv[6]))
massAcc = float(sys.argv[10])
fragThresh= 0.01
offset = .7
useIso = bool(int(sys.argv[7]))
threshold = 0
lam = float(sys.argv[8])
unknowns = bool(int(sys.argv[9]))

if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    decID = DecoID.DecoID(False, libFile, useAuto, numCores, label=str(lam) + "_")

    decID.readData(filename, 2, usePeaks, DDA, massAcc,offset,peakDefinitions=peakFile)
    if unknowns:
        decID.identifyUnknowns()
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold)
