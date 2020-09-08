from DecoID.DecoID import DecoID
libFile = "../databases/mzCloud_reference.db"
key = open("../housekeeping/mzCloud_api_key.txt").readline().rstrip()
useAuto = True
numCores = 15#int(sys.argv[2])
file = "../exampleData/IROA_P1-6_DIA_test_pos1.mzML"#sys.argv[1]

usePeaks = True
DDA = False
massAcc = 10
fragThresh= 0.01
offset = 25
useIso = False
threshold = 0
lam = 0
rtTol = float("inf")
if __name__ == '__main__':
    #multiprocessing.set_start_method("spawn")
    decID = DecoID(libFile, useAuto, numCores,api_key=key)
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions="../exampleData/IROA_p1-6_peak_table_pos_v3.csv")
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)
