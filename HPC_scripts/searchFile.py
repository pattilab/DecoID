import sys
from DecoID import DecoID
import multiprocessing

filename = sys.argv[2]
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
libFile = sys.argv[1]
lab = sys.argv[11]
frag_cutoff = float(sys.argv[12])
rtTol = float(sys.argv[13])

if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    decID = DecoID.DecoID(libFile,"none",numCores,2,label="_"+lab+"_"+str(lam) + "_")
    decID.readData(filename,2,usePeaks,DDA,massAcc,offset,peakDefinitions=peakFile,frag_cutoff=frag_cutoff)
    if unknowns:
        decID.identifyUnknowns(lam,iso=useIso,rtTol=rtTol)
    decID.searchSpectra("y",lam,iso=useIso,threshold=threshold,rtTol=rtTol)
