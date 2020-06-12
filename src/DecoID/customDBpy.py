import os

import sys
if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
elif __file__:
    application_path = os.path.dirname(__file__)


def getCanidateCompoundTrees(mode, upperBound,lowerBound, isotope=False, library="none",cache="none",resolution=2):

    possCompounds = set()
    possIsotopes = set()

    #get compounds in isolation window
    possCompoundsTemp = {(library[mode][x]["cpdID"],library[mode][x]["cpdID"],library[mode][x]["name"],library[mode][x]["id"],x): library[mode][x]["m/z"] for x in
                    library[mode] if  library[mode][x]["m/z"] >= lowerBound and  library[mode][x]["m/z"] <= upperBound}
    possCompoundsTemp = [(str(x[0]), x[1], possCompoundsTemp[x], x[2],x[3],x[4]) for x in possCompoundsTemp]

    possCompounds = possCompounds.union({tuple(x) for x in possCompoundsTemp})

    #get isotopes if necessary
    if isotope:
        possIsotopesTemp = {(library[mode][x]["cpdID"],library[mode][x]["cpdID"],library[mode][x]["name"],library[mode][x]["id"],x): library[mode][x]["m/z"] for x in
                    library[mode] if  library[mode][x]["m/z"] >= lowerBound - 1.00335 and  library[mode][x]["m/z"] <= lowerBound}
        possIsotopesTemp = [(str(x[0]),x[1],possIsotopesTemp[x],x[2],x[3],x[4]) for x in possIsotopesTemp]
        possIsotopes = possIsotopes.union({tuple(x) for x in possIsotopesTemp})

    trees = {x[:-1]:library[mode][x[-1]]["spectrum"] for x in possIsotopes.union(possCompounds)}
    uniqueCPDs = list(set([tuple(x[:-1]) for x in trees]))
    trees2 = {c:{} for c in uniqueCPDs}
    for x in trees:
        trees2[tuple(x[:-1])][x[-1]] = trees[x]
    possCompounds = {key[:-2] for key in possCompounds}
    possIsotopes = {key[:-2] for key in possIsotopes}
    return trees2,possCompounds,possIsotopes

