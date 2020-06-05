import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src/DecoID/"))
from DecoID import DecoID
import multiprocessing

filename = sys.argv[1]#"../exampleData/Asp-40uM_Mal1uM_1Da_NCE30_750ms_NCE50"
filekey = sys.argv[2]
lam = float(sys.argv[4])

if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")
    fragThresh = 0.01
    useIso = False
    threshold = 0

    decIDRec = DecoID.DecoID.fromDill(filename + "_" + filekey + "_.dill")
    decIDRec.numCores = int(sys.argv[3])
    decIDRec.searchSpectra("y",lam, fragThresh, useIso, threshold)


