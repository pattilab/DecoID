import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src/"))
import DecoID
import multiprocessing


useAuto = True
filename = sys.argv[1]

numFiles = int(sys.argv[2])

if __name__ == '__main__':
    multiprocessing.set_start_method("spawn")

    baseName = filename
    endings = ["_" + str(x) + "_" for x in range(numFiles)]
    filenames = [filename + end for end in endings]

    DecoID.DecoID.combineResults(filenames, filename)