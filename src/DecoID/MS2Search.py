import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import scipy.optimize as opt
import sklearn.linear_model as linModel

import os
from pyteomics import mzml
import sys
import sklearn.utils._cython_blas

if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
elif __file__:
    application_path = os.path.dirname(__file__)

from bisect import bisect_left

MAXMASS = 3000


def collapseAsNeeded(spectra,targetRes):
    tempSpec = {np.round(x,targetRes):0 for x in spectra}
    for key,val in spectra.items():
        tempSpec[np.round(key,targetRes)] += val
    return tempSpec

def solveSystem(S,o,resPenalty,maxQuant=1000):
    return deconvolveLASSO(np.transpose(S), [[x] for x in o], [0 for _ in S], [maxQuant for _ in S], resPenalty=resPenalty)

def normalizeSpectra(spectra,method="sum"):
    if method == "sum": maxSpec = np.sum(spectra)
    else: maxSpec = np.max(spectra)
    if np.isinf(np.power(np.float(maxSpec/maxSpec),2)) or np.isnan(np.power(np.float(maxSpec/maxSpec),2)):
        return [0.0 for x in spectra]
    normSpec = [x/maxSpec for x in spectra]
    if np.max(normSpec) > 1.1:
        return [0.0 for x in spectra]
    return [x/maxSpec for x in spectra]

def isPossible(query,candidate,threshold=.75):
    if len(candidate.keys()) == 0:
        return False
    if len(set(query.keys()).intersection(set(candidate.keys())))/float(len(candidate.keys())) > threshold:
        return True
    else:
        return False


def findMaxL1Norm(A,b):
    func = lambda x:np.dot(np.transpose(np.subtract(b,np.dot(A,x*np.ones((len(A[0]),1))))),np.subtract(b,np.dot(A,x*np.ones((len(A[0]),1)))))
    sol = opt.minimize_scalar(func,options={"maxiter":100},tol=1e-1)
    val = sol.x
    return val*len(A[0])

def checkPSD(x):
    delta = .1 * np.mean(flatten(x))
    try:
        np.linalg.cholesky(x)
        return True
    except:
        if np.all(np.linalg.eigvals(x) >= 0):
            return True
        else:
            return False

def deconvolveLASSO(A, b, lb, ub, resPenalty=10,numRepeats = 10):

    scores = []
    sparseQ = True
    if resPenalty == 0:
        sparseQ = False
        resPenalty = 1
    if np.isinf(resPenalty):
        return [[0 for _ in A[0]],0.0]
    if sparseQ:
        model = linModel.Lasso(alpha=resPenalty, fit_intercept=False, positive=True)
        A = np.asfortranarray(A)
        model.fit(A,b,check_input=False)
        params = [[x] for x in model.coef_]
        foundSpectra = np.dot(A,params)
        s2nR = np.sum(flatten(foundSpectra))/np.sum(np.abs(np.subtract(flatten(foundSpectra),flatten(b))))
        return [flatten(params),s2nR]
    else:
        params,_ = opt.nnls(np.array(A),np.array(flatten(b)),maxiter=1e4)
        foundSpectra = np.dot(A, params)
        s2nR = np.sum(flatten(foundSpectra)) / np.sum(
             np.abs(np.subtract(flatten(foundSpectra), flatten(b))))
        return [flatten(params),s2nR]

def dotProductSpectra(foundSpectra,b):

    if type(foundSpectra) == type(dict()):
        mzs = set(foundSpectra.keys()).intersection(set(b.keys()))
        num = np.sum([foundSpectra[x]*b[x] for x in mzs])
        denom = np.sqrt(np.sum([x**2 for x in foundSpectra.values()])*sum([x**2 for x in b.values()]))
        val = num/denom
        return num/denom
    else:
        b = flatten(b)
        foundSpectra = flatten(foundSpectra)
        val = np.sum([x*y for x,y in zip(b,foundSpectra) if x > 1e-4 or y > 1e-4])/np.sqrt(np.sum([x**2 for x in foundSpectra])*(np.sum(
                [x**2 for x in b])))
    if np.isnan(val): val = 0
    return val

def replaceAllCommas(string):
    while "," in string:
        string = string.replace(",","_")
    return string

def splitList(l, n):
    n = int(np.ceil(len(l)/float(n)))
    return list([l[i:i + n] for i in range(0, len(l), n)])

def flatten(l):
    if len(l) > 0 and type(l[0]) == type(l):
        return [item for sublist in l for item in sublist]
    else:
        return l



def scoreDeconvolution(originalSpectra, matrix, res, metIDs, treeIDs, masses, centerMz, DDA, massAcc):

    """Goal is to take deconvolved aspects of the spectra"""
    numComponents = len([x for x in range(len(res)) if res[x] > 0 if type(metIDs[x]) != type(tuple())])
    numComponents += len(set([treeIDs[x] for x in range(len(res)) if res[x] > 0 and type(metIDs[x]) == type(tuple)]))
    resScores = [[0,[0 for _ in originalSpectra],[0 for _ in originalSpectra],"original",masses[x]] for x in range(len(res))]
    if len(matrix) > 0:
        solvedSpectraTotal = np.dot(np.transpose(matrix),[[x] for x in res])
        differences =  np.subtract(flatten(originalSpectra),flatten(solvedSpectraTotal))/numComponents
    else:
        solvedSpectraTotal = [0 for _ in originalSpectra]
        differences = flatten(originalSpectra)



    components = []
    uniqueIsotope = {(metIDs[x][0],t):{"spectrum":-1,"indices":[],"mass":-1,"abundance":0} for x,t in zip(range(len(metIDs)),treeIDs) if type(metIDs[x]) == type(tuple())}
    allUniqueIsotope = {(metIDs[x][0],t):{"spectrum":-1,"indices":[],"mass":-1,"abundance":0} for x,t in zip(range(len(metIDs)),treeIDs) if type(metIDs[x]) == type(tuple())}

    componentsMass = []
    #Form components
    componentName = []
    componentAbudance = []
    allComponents = []
    allComponentMass = []
    allComponentName = []
    allComponentAbundance = []

    resSum = np.sum(res)
    if len([x for x in res if x > 0]) > 1:
        for x in range(len(res)):
            if res[x] > 0:
                if (abs(masses[x]-centerMz)/centerMz*1e6 < massAcc):
                    if type(metIDs[x]) != type(tuple()):
                        components.append([max([res[x] * m + d, 0]) for m, d in zip(matrix[x], differences)])
                        componentsMass.append(masses[x])
                        componentName.append(str(metIDs[x]))
                        componentAbudance.append(res[x])
                    else:
                        uniqueIsotope[(metIDs[x][0],treeIDs[x])]["spectrum"] = [res[x] * m for m in matrix[x]]
                        uniqueIsotope[(metIDs[x][0],treeIDs[x])]["mass"] = masses[x]
                        uniqueIsotope[(metIDs[x][0], treeIDs[x])]["indices"].append(x)
                        uniqueIsotope[(metIDs[x][0], treeIDs[x])]["abundance"] += res[x]


                # closestFrag = [[abs(masses[x]-f),f] for f in fragments]
                # closestFrag.sort(key=lambda x: x[0])
                # closestFrag = closestFrag[0][1]
                if type(metIDs[x]) != type(tuple()):
                    allComponents.append([max([res[x] * m + d, 0]) for m, d in zip(matrix[x], differences)])
                    allComponentMass.append(masses[x])
                    allComponentName.append(str(metIDs[x]))
                    allComponentAbundance.append(res[x])
                else:
                    allUniqueIsotope[(metIDs[x][0], treeIDs[x])]["spectrum"] = [res[x] * m for m in matrix[x]]
                    allUniqueIsotope[(metIDs[x][0],treeIDs[x])]["mass"] = masses[x]
                    allUniqueIsotope[(metIDs[x][0], treeIDs[x])]["indices"].append(x)
                    allUniqueIsotope[(metIDs[x][0], treeIDs[x])]["abundance"] += res[x]



        for x in uniqueIsotope:
            if len(uniqueIsotope[x]["indices"]) > 0:
                components.append([max([p+d,0]) for p,d in zip(uniqueIsotope[x]["spectrum"],differences)])
                componentsMass.append(uniqueIsotope[x]["mass"])
                componentName.append(x[0] + " M+1")
                truePercentage = sum([res[i] for i in uniqueIsotope[x]["indices"]])/len(uniqueIsotope[x]["indices"])
                componentAbudance.append(truePercentage)
                for ind in uniqueIsotope[x]["indices"]:
                    matrix[ind] = uniqueIsotope[x]["spectrum"]
                    res[ind] = truePercentage

        for x in allUniqueIsotope:
            if len(allUniqueIsotope[x]["indices"]) > 0:
                allComponents.append([max([p+d,0]) for p,d in zip(allUniqueIsotope[x]["spectrum"],differences)])
                allComponentMass.append(allUniqueIsotope[x]["mass"])
                allComponentName.append(x[0] + " M+1")
                truePercentage = sum([res[i] for i in allUniqueIsotope[x]["indices"]])/len(allUniqueIsotope[x]["indices"])
                allComponentAbundance.append(truePercentage)
                for ind in allUniqueIsotope[x]["indices"]:
                    matrix[ind] = allUniqueIsotope[x]["spectrum"]
                    res[ind] = truePercentage

    components = {(name,m,ab/resSum):s for m,x,s,name,ab in zip(componentsMass,range(len(componentsMass)),components,componentName,componentAbudance)}
    allComponents = {(name,m,ab/resSum):s for m,x,s,name,ab in zip(allComponentMass,range(len(allComponentMass)),allComponents,allComponentName,allComponentAbundance)}


    #components.append(originalSpectra)
    #componentsMass.append(centerMz)
    components[("original",centerMz,0.0)] = originalSpectra
    #components.append([max([d,0.0]) for d in differences])
    #componentsMass.append(centerMz)
    if len([x for x in res if x > 0]) > 1:
        components[("residual",centerMz,0.0)] = [max([d,0.0]) for d in differences]

    allComponents[("original",centerMz,0.0)] = originalSpectra
    #components.append([max([d,0.0]) for d in differences])
    #componentsMass.append(centerMz)
    if len([x for x in res if x > 0]) > 1:
        allComponents[("residual",centerMz,0.0)] = [max([d,0.0]) for d in differences]

    for x in range(len(res)):
        for comp,mz,ab in components:
            if abs(mz-masses[x])/masses[x]*1e6 < massAcc:
                dp = 100*dotProductSpectra(matrix[x],components[(comp,mz,ab)])
                if dp > resScores[x][0]:
                    resScores[x] = [dp,components[(comp,mz,ab)],matrix[x],comp,centerMz]
    return resScores,allComponents

def safeDivide(num,denom):
    if abs(num) < 1e-12 and abs(denom) < 1e-12:
        return 0.0
    elif abs(num) >= 1e-12 and abs(denom) < 1e-12:
        return np.inf
    else:
        return num/denom



def splitList(l, n):
    n = int(np.ceil(len(l)/float(n)))
    return list([l[i:i + n] for i in range(0, len(l), n)])

def createM1SpectrumfromM0(spectra,massIncrease = 1.003,numFrags=5):

    if len(spectra) > 0:
        newSpectra = []
        frags = [f for f in spectra]
        frags.sort(key=lambda x:spectra[x],reverse=True)
        if len(frags) > numFrags:
            frags = frags[:numFrags]
        for frag in frags:
            tempMzs = dict(spectra)
            tempMzs[frag] = 0
            tempMzs[frag+massIncrease] = spectra[frag]
            newSpectra.append(tempMzs)
        sumSpectra = newSpectra[0]
        for spec in newSpectra[1:]:
            for m in spec:
                if m in sumSpectra:
                    sumSpectra[m] += spec[m]
                else:
                    sumSpectra[m] = spec[m]
    else:
        newSpectra = [{}]
        sumSpectra = {}
    return newSpectra, sumSpectra


def cleanIsotopes(result):
    isotopes = set([(x[0][0],x[1],x[2],x[4]) for x in result if type(x[0]) == type(tuple())])
    update = {(key[3],key[0],key[1],key[2],key[4]):val for key,val in result.items() if type(key[0]) != type(tuple())}
    for iso in isotopes:
        scores = []
        comp = 0.0
        for metID,tree,mass,name,c in result:
            if type(metID) == type(tuple()) and metID[0] == iso[0] and iso[1] == tree:
                scores.append(result[(metID,tree,mass,name,c)][0])
                metOfInterest = metID
                treeOfInterest = tree
                massOfInterest = mass
                nameOfInterest = name
                compOfInterest = c
                comp = c
        scores = np.mean(scores)
        update[(str(nameOfInterest[0]) + " (M+1)",iso[0],iso[1],iso[2],comp)] = [scores] + result[(metOfInterest,treeOfInterest,massOfInterest,nameOfInterest,compOfInterest)][1:]
    return update


def collapse1Fold(spec):
    inc = 10
    ind = 0
    newVect = []
    #print(len(spec))
    while ind+inc < len(spec):
        newVect.append(sum(spec[ind:ind+inc]))
        ind += inc
    newVect.append(sum(newVect[ind:]))
    #print(len(newVect))
    return newVect

def clusterSpectraByMzs(mzs,ids,rts,groups,mzTol = 5):

    prelimClusters = {}
    clusterGroups = {}
    for x in range(len(mzs)):
        good = True
        for x2 in prelimClusters:
            if 1e6*abs(x2 - mzs[x])/mzs[x] <= mzTol or groups[x] in clusterGroups[x2]:
                clusterGroups[x2].append(groups[x])
                prelimClusters[x2].append([mzs[x],ids[x]])
                good = False
                break
        if good:
            prelimClusters[mzs[x]] = [[mzs[x],ids[x]]]
            clusterGroups[mzs[x]] = [groups[x]]

    return prelimClusters

def clusterSpectra(samples,peakData):
    groups = list(set([s["group"] for s in samples]))
    clusters = {g:[] for g in groups}
    for s in samples:
        for g in groups:
            if s["group"] == g:
                clusters[g].append(s)
                break
    clustFormatted = []
    for clust in clusters:
        spec = clusters[clust][0]["spectra"]

        if len(clusters[clust]) > 1:
            for spec2 in clusters[clust][1:]:
                spec = mergeSpectrum(spec,spec2["spectra"])
        clustFormatted.append({"group":clust,"spectrum":spec ,"m/z":np.mean([x["center m/z"] for x in clusters[clust]]),"rtWindow":[peakData.at[clust,"rt_start"],peakData.at[clust,"rt_end"]]})



        #def mergeCluster(cluster1,cluster2):
        #     return {"mergedSpectrum": mergeSpectrum(cluster1["mergedSpectrum"],cluster2["mergedSpectrum"]),
        #             "allSpec":cluster1["allSpec"]+cluster2["allSpec"],"ids":cluster1["ids"]+cluster2["ids"],
        #             "mzs": cluster1["mzs"]+cluster2["mzs"],"rts":cluster1["rts"]+cluster2["rts"]}
        #
        # ppmWindow = lambda m: mzTol * 1e-6 * m
        #
        # uniqueMzs = list(set(mzs))
        # absUniqueMzs = []
        # for x in range(len(uniqueMzs)):
        #     win = ppmWindow(uniqueMzs[x])
        #     good = True
        #     for x2 in absUniqueMzs:
        #         if abs(x2 - uniqueMzs[x]) <= win:
        #             good = False
        #             break
        #     if good:
        #         absUniqueMzs.append(np.round(uniqueMzs[x],6))
        #
        #
        # prelimClusters = {mz:[] for mz in absUniqueMzs}
        #
        # for mz in absUniqueMzs:
        #     win = ppmWindow(mz)
        #     toRemove = []
        #     for m,r,o,i in zip(mzs,rts,order,range(len(order))):
        #         if abs(mz-m) <= win:
        #             prelimClusters[mz].append([o,m,r])
        #             toRemove.append(i)
        #
        #     mzs = [mzs[x] for x in range(len(order)) if x not in toRemove]
        #     rts = [rts[x] for x in range(len(order)) if x not in toRemove]
        #     order = [order[x] for x in range(len(order)) if x not in toRemove]
        # totalClusters = []
        # for c in prelimClusters:
        #     tempClust = [{"mergedSpectrum":spectra[id],"allSpec":[spectra[id]],"ids":[id],"mzs":[mz],"rts":[rt]} for id,mz,rt in prelimClusters[c]]
        #     delta = (1-minScores)/numRounds
        #     thresh = 1.0
        #     numSame = 0
        #     #rd.seed(1300)
        #
        #     for r in range(3*numRounds):
        #         if thresh >= minScores:
        #             thresh = thresh - delta
        #         rd.shuffle(tempClust)
        #         numBefore = len(tempClust)
        #         index = 0
        #         for c in deepcopy(tempClust):
        #             index+=1
        #             for cp in range(tempClust.index(c)):
        #                 if dotProductSpectra(c["mergedSpectrum"][0],tempClust[cp]["mergedSpectrum"][0]) >= .5:
        #                     tempClust[cp] = mergeCluster(tempClust[cp],c)
        #                     tempClust.remove(c)
        #                     break
        #
        #             #clusters = [clusters[x] for x in range(len(clusters)) if x not in toRemove]
        #         if len(tempClust) == numBefore and abs(thresh-minScores) < .01:
        #             numSame += 1
        #             if numSame == 2:
        #                 break
        #         else:
        #             numSame = 0
        #     #print(numRounds)
        #     #rd.seed(datetime.now())
        #     totalClusters += tempClust
    return clustFormatted


def mergeSpectrum(consensous,toAdd):
    def getVal(dict,key):
        if key not in dict:
            return 0.0
        else:
            return dict[key]
    keys = set(list(consensous.keys()) + list(toAdd.keys()))
    new = {x:getVal(consensous,x)+getVal(toAdd,x) for x in keys}
    return new

def findMax(mzs,ms1,rtStart,allRts):
    getInt = lambda rt: np.sum([ms1[rt][x] for x in mzs if x in ms1[rt]])
    rInd = allRts.index(rtStart)
    grad = lambda ind,d: (getInt(allRts[ind+d]) - getInt(allRts[ind]))/(allRts[ind+d]-allRts[ind])#,(getInt(allRts[ind]) - getInt(allRts[ind-1]))/(allRts[ind]-allRts[ind-1])])
    while rInd >= len(allRts) -1:
        rInd -= 1
    while rInd <= 1:
        rInd += 1
    pk = grad(rInd,1)
    pkOld = pk
    while(pk * pkOld > 0 and rInd-1 > 1 and rInd+1 < len(allRts)-1):
        pkOld = pk
        if pk > 0:
            rInd += 1
            try:
                pk = grad(rInd,1)
            except:
                print("1",rInd,len(allRts))
        else:
            try:
                rInd -= 1
                pk = grad(rInd,-1)
            except:
                print("-1",rInd,len(allRts))
        #pk = grad(rInd)
        #plt.plot([allRts[rInd],allRts[rInd]],[0,1e7])

    if pk * pkOld <= 0:
        if pkOld > 0:
            rInd -= 1
        else:
            rInd += 1


    return getInt(allRts[rInd]),allRts[rInd]




def extractFeatures(mzs,rts,ids,ms1,mError = 5):
    allRt = list(ms1.keys())
    allRt.sort()
    chromatograms = {}
    allMzList = list()
    for s in ms1:
        allMzList += [x for x in ms1[s] if ms1[s][x] > 0]
    allMzList = np.array(list(set(allMzList)))

    #print("found mzs")
    for id,mz,rt in zip(ids,mzs,rts):
        #plt.figure()
        temp = [[r,abs(r-rt)] for r in allRt]
        temp.sort(key=lambda x: x[1])
        rtRep = temp[0][0]
        mzOI = np.extract(np.abs((10**6)*(allMzList-mz))/mz < mError,allMzList)
        #print("foundmzOI")
        getInt = lambda rt: np.sum([ms1[rt][x] for x in mzOI if x in ms1[rt]])
        maxInt,apex = findMax(mzOI,ms1,rtRep,allRt)
        # maxInt = np.max(list(tempMs1.values()))
        #print("XCR Done")
        r = allRt.index(rtRep)
        tempInts = getInt(allRt[r])
        while(tempInts > .001*maxInt and r > 0 and r < len(allRt)-1):
            r -= 1
            tempInts = getInt(allRt[r])

        rtMin = allRt[r]
        #print("rt min done")
        r = allRt.index(rtRep)
        tempInts = getInt(allRt[r])
        while(tempInts > .001*maxInt and r > 0 and r < len(allRt)-1):
            r += 1
            tempInts = getInt(allRt[r])


        rtMax = allRt[r]
        #print("rtMaxDone")
        chromatograms[id] = [rtMin,rtMax]
        rtMinBound = rtMin-.5*(rtMax-rtMin)
        rtMaxBound = rtMax+.5*(rtMax-rtMin)
        #tempMs1 = [[key,getInt(key)] for key in allRt if key > rtMinBound and key < rtMaxBound]


        # plt.figure()
        # plt.plot([x[0] for x in tempMs1],[x[1] for x in tempMs1])
        # plt.plot([rt,rt],[0,maxInt])
        # plt.plot([rtMin,rtMin],[0,maxInt])
        # plt.plot([rtMax,rtMax],[0,maxInt])
        # plt.plot([apex,apex],[0,maxInt])
        # plt.xlim([rtMinBound,rtMaxBound])

        # plt.figure()
        # plt.plot([allRt[x] for x in range(len(allRt)-1) if allRt[x] > rtMinBound and allRt[x] < rtMaxBound],
        #          [(tempMs1[allRt[x+1]] - tempMs1[allRt[x]])/(allRt[x+1] - allRt[x]) for x in range(len(allRt)) if allRt[x] > rtMinBound and allRt[x] < rtMaxBound])
        # plt.plot([rt, rt], [0, maxInt])
        # plt.plot([rtMin, rtMin], [0, maxInt])
        # plt.plot([rtMax, rtMax], [0, maxInt])
        # plt.xlim([rtMinBound, rtMaxBound])

        #print("plots done")
    #plt.show()
    return chromatograms


def getMatricesForGroup(trees,spectra,possCompounds,possIsotopes,isotope,res,clusters):

    spectrum = dict(spectra[0])
    for spec in spectra[1:]:
        for m,i in spec.items():
            if m in spectrum:
                spectrum[m] += i
            else:
                spectrum[m] = i

    compoundDict = pullMostSimilarSpectra({key:val for key,val in trees.items() if key in possCompounds},spectrum)
    compoundDict = {key: [val[0], collapseAsNeeded(val[1],res)] for key,val in compoundDict.items()}#convertSpectra2Vector(val[0], MAXMASS, res)] for key, val in
                    #compoundDictFull.items()}

    keys = list(compoundDict.keys())

    if isotope:
        isotopeDictFull = pullMostSimilarSpectraIsotopologues({key:val for key,val in trees.items() if key in possIsotopes}, spectrum)

        isotopeDict2 = {}
        for key, val in isotopeDictFull.items():
            isotopeSpectras = val[1]
            for spec, num in zip(isotopeSpectras, range(len(isotopeSpectras))):
                tempKey = ((str(key[0]), num), key[1], key[2] + 1.003, (key[3], num))
                isotopeDict2[tempKey] = [key[1], spec]
        isotopeKeys = list(isotopeDict2.keys())
        isotopeDict = {key: [val[0], collapseAsNeeded(val[1], res)] for key, val in isotopeDict2.items()}


    if type(clusters) != type(None):
        clusterKeys = list(clusters.keys())

    metIDs = [x[0] for x in keys]
    spectraTrees = [x[1] for x in keys]
    masses = [x[2] for x in keys]
    metNames = [x[3] for x in keys]
    spectraIDs = [compoundDict[x][0] for x in keys]
    matrix = [compoundDict[x][1] for x in keys]

    if isotope:
        metIDs += [x[0] for x in isotopeKeys]
        masses += [x[2] for x in isotopeKeys]
        spectraTrees += [x[1] for x in isotopeKeys]
        spectraIDs += [isotopeDict[x][0] for x in isotopeKeys]
        metNames += [x[3] for x in isotopeKeys]
        matrix += [isotopeDict[x][1] for x in isotopeKeys]

    if type(clusters) != type(None):
        metIDs += [str(key[0]) + "|" + str(np.mean(key[2])) for key in clusterKeys]
        masses += [key[0] for key in clusterKeys]
        spectraTrees += [-1 for _ in clusterKeys]
        spectraIDs += [-1 for _ in clusterKeys]
        metNames += [
            "Predicted Analyte: m/z = " + str(np.round(key[0], 5)) + " rt_Range: [" + str(key[1]) + "-" + str(
                key[2]) + "]" for key in clusterKeys]
        matrix += [clusters[key] for key in clusterKeys]

        # convert dicts to lists for deconvolution

    coeffs = []
    for m in metIDs:
        if type(m) != type(tuple()):
            coeffs.append(1.0)
        else:
            coeffs.append(1.0)
            # coeffs.append(len([m1 for m1 in mets if type(m1) == type(tuple()) and m1[0] == m[0]]))
    indexRef = [np.round(x * 10 ** (-1 * res), res) for x in list(range(int(MAXMASS * 10 ** res)))]
    indices = list(set(flatten([list(spectrum.keys())] + [list(m.keys()) for m in matrix])))
    indices.sort()

    indicesAll = [indexRef.index(np.round(x, res)) for x in indices]
    matrix = [[s * c for s in normalizeSpectra([getVal(m, x) for x in indices])] for m, c in
                      zip(matrix, coeffs)]
    reduceSpec = [[getVal(spec, x) for x in indices] for spec in spectra]

    return metIDs, spectraTrees, spectraIDs, matrix, masses, metNames, indicesAll, reduceSpec

def pullMostSimilarSpectraIsotopologues(trees,spectra):
    returnDict = {}
    for tree,ms2Scans in trees.items():
        if len(ms2Scans) > 0:
            generatedIsotopologus = {id:createM1SpectrumfromM0(ms2) for id,ms2 in ms2Scans.items()}
            temp = [[id,ms2[0],dotProductSpectra(spectra,ms2[1])] for id,ms2 in generatedIsotopologus.items()]
            temp.sort(key=lambda x:x[2],reverse=True)
            returnDict[tree] = temp[0][:2]
    return returnDict


def pullMostSimilarSpectra(trees,spectra):
    returnDict = {}
    for tree,ms2Scans in trees.items():
        if len(ms2Scans) > 0:
            temp = [[id,ms2,dotProductSpectra(spectra,ms2)] for id,ms2 in ms2Scans.items()]
            temp.sort(key=lambda x:x[2],reverse=True)
            returnDict[tree] = temp[0][:2]
    return returnDict


# get value from dictionary or return 0 if it is not in the dictionary
def getVal(dict,key):
    if key in dict:
        return dict[key]
    else:
        return 0.0


# determine if m/z is within ppmThresh of an m/z in fragments
def inScan(fragments,candidate_mz,ppmThresh=10):
    if type(fragments) != type(str()):
        try:
            if any(abs(candidate_mz-x)/candidate_mz/(1e-6) < ppmThresh for x in fragments):
                return True
            else:
                return False
        except:
            print(candidate_mz,fragments)
    else:
        return True

def inScanIso(fragments,candidate_mz,ppmThresh=10):
    if type(fragments) != type(str()):
        if any(abs(candidate_mz - x) / candidate_mz / (1e-6) < ppmThresh for x in fragments) and any(abs(candidate_mz-1.003 - x)/(candidate_mz-1.003)/(1e-6) < ppmThresh for x in fragments):
            return True
        else:
            return False
    else:
        return True


#https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before


def readRawDataFile(filename, maxMass, resolution, useMS1, ppmWidth = 50,offset=0.65,tic_cutoff=5):
   try:
        delete = False
        if ".mzML" not in filename:
            toRemove = filename.split("/")[-1]
            command = "msconvert " + filename + " --outdir " + filename.replace(toRemove,"") + ' --filter "peakPicking true 1-" --ignoreUnknownInstrumentError > junk.txt'
            #command = "msconvert " + filename + " --outdir " + filename.replace(toRemove,"") + ' --ignoreUnknownInstrumentError > junk.txt'
            os.system(command)
            filename = filename.replace(".raw", "")
            filename = filename.replace('"', "")
            filename = filename+".mzML"
            delete = True
        reader = mzml.read(filename.replace('"', ""))
        result = []
        ms1Scans = {}
        scanIDCount = 0
        for temp in reader:
            if temp['ms level'] == 2:
                    tic = np.log10(float(temp["total ion current"]))
                    if tic >= tic_cutoff:
                        id = scanIDCount#temp["id"].split()[-1].split("=")[-1]
                        scanIDCount += 1
                        #centerMz = temp["precursorList"]["precursor"][0]["isolationWindow"]['isolation window target m/z']
                        centerMz = temp["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
                        try:
                            lowerBound = centerMz - temp["precursorList"]["precursor"][0]["isolationWindow"]['isolation window lower offset']
                            upperBound = centerMz + temp["precursorList"]["precursor"][0]["isolationWindow"]['isolation window upper offset']
                        except:
                            lowerBound = centerMz - offset
                            upperBound = centerMz + offset
                        #filter = temp["scanList"]["scan"][0]["filter string"]
                        #acquisitionMode = filter.split("@")[0].split()[1]
                        if 'positive scan' in temp:
                            acquisitionMode = "Positive"
                        else:
                            acquisitionMode = "Negative"
                        #settings = filter.split("@")[1].split()[0]
                        #fragmentMode = settings[:3]
                        #NCE = float(settings[3:])
                        rt = temp["scanList"]["scan"][0]["scan start time"]
                        mzs = list(zip(temp["m/z array"],temp["intensity array"]))
                        tempSpecs = []

                        spectra ={np.round(x[0],resolution):0 for x in mzs}
                        for x,y in mzs: spectra[np.round(x,resolution)] += y

                        result.append({"id":id,"spectra":spectra,"mode":acquisitionMode,"center m/z":
                                           centerMz,"lower m/z":lowerBound,"higher m/z":upperBound,"rt":rt,"signal":tic})
            elif useMS1 and temp['ms level'] == 1:
                ms1Scans[temp["scanList"]["scan"][0]["scan start time"]] = {mz: i for mz, i in zip(temp["m/z array"], temp["intensity array"])}
        reader.close()
        if len(ms1Scans) > 0 and useMS1:
            rts = list(ms1Scans.keys())
            rts.sort()
            for samp in range(len(result)):
                isoWidth = result[samp]["center m/z"] - result[samp]["lower m/z"]

                mzScan = ms1Scans[takeClosest(rts,result[samp]["rt"])]

                peaksOfInterest = [[m, i] for m, i in mzScan.items() if m >= result[samp]["center m/z"] - isoWidth - 1.3
                                   and m <= result[samp]["center m/z"] + isoWidth]
                peaksOfInterest.sort(key=lambda x: x[0])
                precursorPeaks = [x for x in peaksOfInterest if
                                  abs(x[0] - result[samp]["center m/z"]) * 1e6 / result[samp]["center m/z"] <= ppmWidth and x[1] > 0]
                chimericPeaks = [x for x in peaksOfInterest if
                                 abs(x[0] - result[samp]["center m/z"]) * 1e6 / result[samp]["center m/z"] > ppmWidth and
                                 x[0] >= result[samp]["center m/z"] - isoWidth and x[1] > 0]
                if len(chimericPeaks) > 0 or len(precursorPeaks) > 0:
                    # result[samp]["percentContamination"] = min([1, max([0, simps([x[1] for x in chimericPeaks],
                    #                                                               [x[0] for x in chimericPeaks],
                    #                                                               even="avg") / (
                    #                                                              simps([x[1] for x in chimericPeaks],
                    #                                                                    [x[0] for x in chimericPeaks],
                    #                                                                    even="avg") + simps(
                    #                                                          [x[1] for x in precursorPeaks],
                    #                                                          [x[0] for x in precursorPeaks],
                    #                                                          even="avg"))]), ])
                    result[samp]["percentContamination"] = np.sum([x[1] for x in chimericPeaks])/(np.sum([x[1] for x in chimericPeaks]) + np.sum([x[1] for x in precursorPeaks]))
                # elif len(chimericPeaks) > 0:
                #     result[samp]["percentContamination"] = 1.0
                else:
                    result[samp]["percentContamination"] = 0.0
                result[samp]["fragments"] = [x[0] for x in peaksOfInterest if x[1] > 1e-6]
                result[samp]["ms1"] = [x for x in peaksOfInterest if x[0] >= result[samp]["center m/z"] - isoWidth]

            if delete:
                #os.remove(filename)
                os.remove("junk.txt")
        return result,ms1Scans

   except:
        print(sys.exc_info())
        print(filename + " does not exist or is ill-formatted")
        return -1,-1


def createDictFromString(string, resolution):
    specDict = {}
    for val in string.split():
        mz = np.round(float(val.split(":")[0]), resolution)
        if mz in specDict:
            specDict[mz] += float(val.split(":")[1])
        else:
            specDict[mz] = float(val.split(":")[1])
    return specDict