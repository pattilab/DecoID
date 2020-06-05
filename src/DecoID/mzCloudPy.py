import grequests
import json
import os
from keys import *
import pickle as pkl

import sys
if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
    MZCOMPOUNDTREELINK = {"reference": pkl.load(
        open(os.path.join(application_path, "mzCloudCompound2TreeLinkagereference.pkl"), "rb")),
                          "autoprocessing": pkl.load(open(os.path.join(application_path,
                                                                       "mzCloudCompound2TreeLinkageautoprocessing.pkl"),
                                                          "rb"))}
elif __file__:
    application_path = os.path.dirname(__file__)
    MZCOMPOUNDTREELINK = {"reference": pkl.load(
        open(os.path.join(application_path, "../../databases/mzCloudCompound2TreeLinkagereference.pkl"), "rb")),
                          "autoprocessing": pkl.load(open(os.path.join(application_path,
                                                                       "../../databases/mzCloudCompound2TreeLinkageautoprocessing.pkl"),
                                                          "rb"))}

CACHEFILE = os.path.join(application_path, 'cacheReq')
from MS2Search import *



timeout = 30
CONCURRENTREQUESTMAX = 1

MAXMASS = 5000



def getCanidateCompoundTrees(mode, upperBound,lowerBound, isotope=False, library=["reference"],cacheReq = "none"):

    #make list of isolation window range
    centerMz = [lowerBound,upperBound]
    centerMz[0] = int(np.floor(centerMz[0]))
    centerMz[1] = int(np.ceil(centerMz[1]))
    possCompoundsTrueAll = set()
    possIsotopesTrueAll = set()
    possCompounds = set()
    possIsotopes = set()
    for lib in library:
        #get compounds in isolation window
        possCompoundsTemp = {tuple(x[:2]) + (x[3],): x[2] for x in
                        flatten([MZCOMPOUNDTREELINK[lib][mode][x] for x in range(centerMz[0] - 1, centerMz[1] + 1)])}
        possCompoundsTemp = [(lib[0] + str(x[0]), x[1], possCompoundsTemp[x], x[2]) for x in possCompoundsTemp if possCompoundsTemp[x] >= lowerBound and possCompoundsTemp[x] <= upperBound]

        possCompounds = possCompounds.union({tuple(x) for x in possCompoundsTemp})

        #get isotopes if necessary
        if isotope:
            possIsotopesTemp = {tuple(x[:2]) + (x[3],):x[2] for x in flatten([MZCOMPOUNDTREELINK[lib][mode][x] for x in range(centerMz[0]-1, centerMz[1])])}
            possIsotopesTemp = [(lib[0] + str(x[0]),x[1],possIsotopesTemp[x],x[2]) for x in possIsotopesTemp if possIsotopesTemp[x] >= lowerBound - 1.00335  and possIsotopesTemp[x] <= lowerBound]

            possIsotopes = possIsotopes.union({tuple(x) for x in possIsotopesTemp})# if x[2] >= centerMzOrig - isolationWidth - 1 and x[2] <= centerMzOrig - isolationWidth}


    trees = getTreesAsync(possIsotopes.union(possCompounds),cache=cacheReq)

    return trees,possCompounds,possIsotopes


"""
take spectra from m/z cloud after being read in by json and convert to dictionary
"""
def convertSpectra2Vector(spectra,maxMass,resolution):
    mzs = {np.round(x["MZ"],resolution):x["Abundance"] for x in spectra["Peaks"]}
    return mzs

def getTreesAsync(trees, calibration="recalibrated",maxRetries = 3,cache="none",maxMass = MAXMASS,resolution = 2):
    libraryDict = {"a":"autoprocessing","r":"reference"}
    requestTuple = []
    keyList = []
    alreadyCached = {}
    toCache = []
    for tree in trees:
        url = SPECTRAURL.replace("TREENUMBER", str(tree[1]))
        library = libraryDict[tree[0][0]]
        url = url.replace("LIBRARY",library)
        querystring = {"stage": "2", "processing": calibration, "peaks": "true"}
        payload = ""
        if type(cache) != type("") and (library,tree) in cache:
            alreadyCached[tree] = cache[(library,tree)]
        else:
            toCache.append([library,tree])
            keyList.append(tree)
            requestTuple.append(grequests.get(url, data=payload, headers=HEADERSSPEC, params=querystring, timeout=timeout))
    requestTuple = tuple(requestTuple)
    responses = grequests.map(requestTuple,size=CONCURRENTREQUESTMAX)
    errors = [key for key in range(len(keyList)) if
              type(responses[key]) == type(None) or "An error has occurred" in responses[key].text or "Service Unavailable" in responses[key].text]
    numTries = 1
    while(len(errors) > 0) and maxRetries > numTries:
        requestTuple = []
        #print("error", len(errors),len(trees))
        for tree in errors:
            url = SPECTRAURL.replace("TREENUMBER", str(keyList[tree][1]))
            url = url.replace("LIBRARY", library)
            querystring = {"stage": "2", "processing": calibration, "peaks": "true"}
            payload = ""
            requestTuple.append(
                grequests.get(url, data=payload, headers=HEADERSSPEC, params=querystring, timeout=timeout))
        requestTuple = tuple(requestTuple)
        responsesNew = grequests.map(requestTuple,size=CONCURRENTREQUESTMAX)
        for response,tree in zip(responsesNew,errors):
            responses[tree] = response
        errors = [key for key in range(len(keyList)) if type(responses[key]) == type(None)or "An error has occurred" in responses[key].text or "Service Unavailable" in responses[key].text]
        numTries+=1
    if numTries >= maxRetries:
        print("Number of retries exceeded to contact m/zCloud. Check network")
    responses = [responses[x] for x in range(len(responses)) if x not in errors]
    keyList = [keyList[x] for x in range(len(keyList)) if x not in errors]
    output = {key:getAllSpectraInTree(reformatSpectraDictList(json.loads(val.text)),maxMass,resolution) for key,val in zip(keyList,responses)}
    #output = {key:reformatSpectraDictList(json.loads(val.text)) for key,val in zip(keyList,responses)}
    if type(cache) != type(""):
        tempDict = {}
        for x in toCache:
            try:
                tempDict[(x[0],x[1])] = output[x[1]]
            except:
                pass
        output.update(alreadyCached)
        chunkSize = 500
        allKeys = list(tempDict.keys())
        if len(allKeys) > 0:
            for smallChunks in splitList(allKeys, max([1, int(len(allKeys) / chunkSize)])):
                cache.update({key: tempDict[key] for key in smallChunks})
        #print(len(tempDict),len(output))
    return output

def reformatSpectraDictList(dictList):
    dat = {x["Id"]:{key:val for key,val in x.items() if key != "Id"} for x in dictList}
    return dat


def getAllSpectraInTree(dat,maxMass,resolution):
    return {id:convertSpectra2Vector(dat[id],maxMass,resolution) for id in dat}


def getCompoundList(page,pageSize = 100,library="reference"):
    url = COMPOUNDURL
    url = url.replace("LIBRARY",library)
    querystring = {"page": str(page), "pageSize": str(pageSize), "newer": "2000-09-03"}
    payload = ""

    response = 1
    while(response==1):
        try:
            response = (grequests.get(url, data=payload, headers=HEADERSCOMP, params=querystring,timeout=timeout),)
            response = grequests.map(response)[0]
            if type(response) == type(None) or "An error has occurred" in response.text or "Service Unavailable" in response.text:
                response = 1
            else:
                dat = reformatSpectraDictList(json.loads(response.text)["Items"])
        except:
            pass
    return dat



def getListofSpectrainTree(trees,calibration="recalibrated",library="reference"):
    response = 1
    requestTuple = []
    treeOrder = list(trees)
    for tree in trees:
        url = SPECTRAURL.replace("TREENUMBER", str(tree[0]))
        url = url.replace("LIBRARY",library)
        querystring = {"stage": "2", "processing": calibration, "peaks": "true"}
        payload = ""
        requestTuple.append(grequests.get(url, data=payload, headers=HEADERSSPEC, params=querystring, timeout=40))
    requestTuple = tuple(requestTuple)
    responses = grequests.map(requestTuple)

    errors = [key for key in range(len(responses)) if
              type(responses[key]) == type(None) or "An error has occurred" in responses[key].text]
    while (len(errors) > 0):
        requestTuple = []
        for tree in [treeOrder[x] for x in errors]:
            url = SPECTRAURL.replace("TREENUMBER", str(tree[0]))
            url = url.replace("LIBRARY", library)
            querystring = {"stage": "2", "processing": calibration, "peaks": "true"}
            payload = ""
            requestTuple.append(grequests.get(url, data=payload, headers=HEADERSSPEC, params=querystring, timeout=40))
        requestTuple = tuple(requestTuple)
        responsesNew = grequests.map(requestTuple)
        for error,res in zip(errors,responsesNew):
            responses[error] = res

        errors = [key for key in range(len(responses)) if
                  type(responses[key]) == type(None) or "An error has occurred" in responses[key].text or "Service Unavailable" in responses[key].text or "unavailable" in responses[key].text]

    output = {tree:reformatSpectraDictList(json.loads(response.text)) for response,tree in zip(responses,treeOrder)}
    return output


def generateCompoundID2SpectralIDIndexedByM_ZStrict(numPerPage = 100,library="reference"):
    #get number of items

    totalCompounds = json.loads(grequests.map((grequests.get(COMPOUNDURL.replace("LIBRARY",library), data="", headers=HEADERSCOMP, params={"page": str(1),
          "pageSize": "5", "newer": "2000-09-03"}),))[0].text)["Total"]
    numPages = int(np.ceil(float(totalCompounds)/numPerPage))
    pos = [[] for _ in range(5000)]
    neg = [[] for _ in range(5000)]
    linkage = {"Positive":pos,"Negative":neg}
    for page in range(1,numPages+1):
        compounds = getCompoundList(page,numPerPage,library=library)
        #print(compounds)
        treeDict = {}
        for comp in compounds:
            treeDict.update({(key,comp,compounds[comp]["SearchCompoundName"]):val for key,val in reformatSpectraDictList(compounds[comp]["SpectralTrees"]).items()})

        allSpectra = getListofSpectrainTree(list(treeDict.keys()))

        for tree in treeDict:
            posValsRounded = []
            posValsTrue = []
            spectras = allSpectra[tree]
            for spectra in spectras:
                rounded = int(np.floor(spectras[spectra]["IsolationWidth"]["XPos"]))
                posValsRounded.append(rounded)
                posValsTrue.append(np.round(spectras[spectra]["IsolationWidth"]["XPos"],7))
            try:
                [linkage[treeDict[tree]["Polarity"]][x].append((tree[1],tree[0],y,tree[2])) for x,y in zip(posValsRounded,posValsTrue)]
            except:
                print(tree,posValsRounded)
        print(float(page)*numPerPage/totalCompounds)
    pkl.dump(linkage,open("mzCloudCompound2TreeLinkage"+library+".pkl","wb"),pkl.HIGHEST_PROTOCOL)