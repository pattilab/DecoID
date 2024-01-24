
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
import sys
import os
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)
from multiprocessing import Queue, Process,Manager,Value,Lock,Pool
from threading import Thread
import time
import gzip
import pickle as pkl
import dill
import pandas as pd
import itertools
import uuid
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from copy import deepcopy
import requests
import json
import scipy.stats as stats
import logging
import os
import IsoSpecPy
import zipfile

if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
    # MZCOMPOUNDTREELINK = {"reference": pkl.load(
    #     open(os.path.join(application_path, "mzCloudCompound2TreeLinkagereference.pkl"), "rb")),
    #     "autoprocessing": pkl.load(open(os.path.join(application_path,
    #                                                  "mzCloudCompound2TreeLinkageautoprocessing.pkl"),
    #                                     "rb"))}
elif __file__:

    application_path = os.path.dirname(__file__)

MZCOMPOUNDTREELINK = {"reference": pkl.load(
    open(os.path.join(application_path, "mzCloudCompound2TreeLinkagereference.pkl"), "rb")),
    "autoprocessing": pkl.load(open(os.path.join(application_path,
                                                 "mzCloudCompound2TreeLinkageautoprocessing.pkl"),
                                    "rb"))}
#uniqueLosses = pkl.load(open(os.path.join(application_path,
#                                                     "uniqueLosses.pkl"),
#                                        "rb"))

path = os.path.join(application_path,"model")

model = -1

logReg = {x:pkl.load(open(os.path.join(application_path,"frac_shared_frag_log_reg_" + x + ".pkl"),"rb")) for x in ["Positive","Negative"]}

timeout = 30
CONCURRENTREQUESTMAX = 1
MAXMASS = 5000
maxMzForPrediction = 500

import molmass
import scipy.optimize as opt
import sklearn.linear_model as linModel

from pyteomics import mzml

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

def normalizeSpectra(spectra, method="sum"):
    if method == "sum": 
        maxSpec = np.sum(spectra)
    else: 
        maxSpec = np.max(spectra)
        
    # Check if maxSpec is zero or NaN
    if maxSpec == 0 or np.isnan(maxSpec):
        return [0.0 for x in spectra]
    
    normSpec = [x/maxSpec for x in spectra]
    if np.max(normSpec) > 1.1:
        return [0.0 for x in spectra]
    return normSpec

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


def dotProductSpectra(foundSpectra,b,mz1=-1,mz2=-1,polarity=-1):
    """
    Computes the normalized dot product (cosine) similarity between two spectra.

    :param foundSpectra: dict or array like, this is the first spectrum
    :param b: dict or array like, this is the second spectrum
    :param mz1: not used for this method
    :param mz2: not used for this method
    :param polarity: not used for this method
    :return: Cosine similarity of the two spectrum in a scale of 0-1
    """
    #if input spectra are dictionaries

    if type(foundSpectra) == type(dict()):
        mzs = set(foundSpectra.keys()).intersection(set(b.keys())) #get shared mzs
        num = np.sum([foundSpectra[x]*b[x] for x in mzs]) #compute num
        denom = np.sqrt(np.sum([x**2 for x in foundSpectra.values()])*sum([x**2 for x in b.values()])) #compute denom
    #if input spectra are lists
    else:
        b = flatten(b) #flatten spectra
        foundSpectra = flatten(foundSpectra)
        num = np.sum([x*y for x,y in zip(b,foundSpectra) if x > 1e-4 or y > 1e-4])
        denom = np.sqrt(np.sum([x**2 for x in foundSpectra])*(np.sum([x**2 for x in b])))

    if denom == 0:
        val = 0
    else:
        val = num/denom
        if np.isnan(val):
            val = 0

    return val

def safeNormalize(vec):
    s = np.sum(vec)
    if s > 1e-5:
        return vec / s
    else:
        return vec

#add shared framgnets
def NNScoring(found,ref,mz1,mz2,polarity):

    #update
    maxMass = 1000
    embeddedSpectra = [np.zeros((maxMass + 1)) for _ in range(2)]
    ind = 0
    for spectrum,_ in zip([found,ref],[mz1,mz2]):
        f = collapseAsNeeded(spectrum,0)
        for mzF, i in f.items():
            mz2_corr = int(np.round(mzF))
            if mz2_corr > maxMass:
                mz2_corr = 0
            embeddedSpectra[ind][mz2_corr] += 1
        ind += 1
    embeddedSpectra = [safeNormalize(x) for x in embeddedSpectra]

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # FATAL
    logging.getLogger('tensorflow').setLevel(logging.FATAL)
    import keras
    if not os.path.isdir(path):
        zipfile.ZipFile(os.path.join(application_path, "model.zip"), "r").extractall(
            os.path.join(application_path, "model"))

    modelName = {"Positive": os.path.join(path, "model_pos"),
             "Negative": os.path.join(path, "model_neg")}

    model = keras.models.load_model(modelName[polarity])

    score = model.predict(np.array([np.concatenate(embeddedSpectra)]))[0][0]

    return score

def sharedFragmentScoring(found,ref,mz1,mz2,polarity):
    if len(found) + len(ref) > 0:
        foundKeys = set({k for k,i in found.items() if i > 1e-5})
        refKeys = set({k for k,i in ref.items() if i > 1e-5})
        val =len(foundKeys.intersection(refKeys))/len(foundKeys.union(refKeys))
    else:
        val = 0
    return logReg[polarity].predict_proba([[val]])[0][1]



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



def scoreDeconvolution(originalSpectra, matrix, res, metIDs, masses, centerMz, rts, massAcc,rtTol,rt,redundancyCheckThresh,scoringFunc,indices,resolution,polarity):

    """Goal is to take deconvolved aspects of the spectra"""
    numComponents = len([x for x in range(len(res)) if res[x] > 0])
    resScores = [[0,[0 for _ in originalSpectra],[0 for _ in originalSpectra],"original",masses[x],rts[x],False] for x in range(len(res))]
    if len(matrix) > 0:
        solvedSpectraTotal = np.dot(np.transpose(matrix),res)
        differences =  np.subtract(flatten(originalSpectra),flatten(solvedSpectraTotal))/numComponents
    else:
        differences = flatten(originalSpectra)

    #Form components
    components = []
    componentsMass = []
    componentRts = []
    componentName = []
    componentAbudance = []

    allComponents = []
    allComponentsMass = []
    allComponentsRts = []
    allComponentsName = []
    allComponentsAbundance = []

    resSum = np.sum(res)
    if len([x for x in res if x > 0]) > 1:
        for x in range(len(res)):
            if res[x] > 0:
                if (abs(masses[x]-centerMz)/centerMz*1e6 < massAcc) and (rts[x] < 0 or abs(rts[x] - rt) < rtTol):
                    components.append([max([res[x] * m + d, 0]) for m, d in zip(matrix[x], differences)])
                    componentsMass.append(masses[x])
                    componentName.append(str(metIDs[x]))
                    componentAbudance.append(res[x]/resSum)
                    componentRts.append(rts[x])

                allComponents.append([max([res[x] * m + d, 0]) for m, d in zip(matrix[x], differences)])
                allComponentsMass.append(masses[x])
                allComponentsName.append(str(metIDs[x]))
                allComponentsAbundance.append(res[x]/resSum)
                allComponentsRts.append(rts[x])


    components = {(name,m,r,ab):s for m,s,name,ab,r in zip(componentsMass,components,componentName,componentAbudance,componentRts)}
    allComponents = {(name,m,r,ab):s for m,s,name,ab,r in zip(allComponentsMass,allComponents,allComponentsName,allComponentsAbundance,allComponentsRts)}

    components[("original",centerMz,rt,0.0)] = originalSpectra
    if len([x for x in res if x > 0]) > 1:
        components[("residual",centerMz,rt,0.0)] = [max([d,0.0]) for d in differences]

    allComponents[("original",centerMz,rt,0.0)] = originalSpectra
    if len([x for x in res if x > 0]) > 1:
        allComponents[("residual",centerMz,rt,0.0)] = [max([d,0.0]) for d in differences]

    redundantComponent = []
    for comp,mz,r,ab in components:
        componentScores = []
        componentThresh = np.inf
        for x in range(len(res)):
            if abs(centerMz-masses[x])/masses[x]*1e6 < massAcc and (rts[x] < 0 or abs(rts[x]-rt) < rtTol):
                #dp = 100*dotProductSpectra(matrix[x],components[(comp,mz,r,ab)])


                spec1 = {}
                spec2 = {}
                for index,i1,i2 in zip(range(len(matrix[x])),matrix[x],components[(comp,mz,r,ab)]):
                    if i1 + i2 > 1e-5:
                        mm = np.round(indices[index] * 10 ** (-1 * resolution), resolution)
                        spec1[mm] = i1
                        spec2[mm] = i2

                dp = 100*scoringFunc(spec1,spec2,masses[x],mz,polarity)

                if dp > resScores[x][0]:
                    resScores[x] = [dp,components[(comp,mz,r,ab)],matrix[x],comp,centerMz,rt,False]
                if metIDs[x] == comp:
                    componentThresh = redundancyCheckThresh * dp
                componentScores.append(dp)
        if len([x for x in componentScores if x > componentThresh]) > 1:
            redundantComponent.append(comp)

    for x in range(len(resScores)):
        if resScores[x][3] in redundantComponent:
            resScores[x][-1] = True




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


def scoreIsotopePattern(ms1,formula,polarity,ppm=10):
    try:
        f1 = molmass.Formula(formula)

        if polarity == "Positive":
            f1 = f1 + molmass.Formula("H")
        else:
            f1 = f1 - molmass.Formula("H")

        formula = f1.formula

        isotopeSpectrum = {m:i for m,i in IsoSpecPy.IsoTotalProb(.99,formula)}

        expected = []
        observed = []
        notFound = list(isotopeSpectrum.keys())

        for mz,i in ms1.items():
            toRemove = []
            found = False
            tmp = 0
            for mz2 in range(len(notFound)):
                if 1e6 * np.abs(notFound[mz2]-mz)/mz < ppm:
                    toRemove.append(mz2)
                    tmp += isotopeSpectrum[notFound[mz2]]
                    found = True
            if found:
                observed.append(i)
                expected.append(tmp)
            notFound = [notFound[x] for x in range(len(notFound)) if x not in toRemove]

        #tmp = 0
        #for mz,i in isotopeSpectrum.items():
        #    for mz2,_ in ms1.items():
        #        if 1e6 * np.abs(mz2-mz)/mz < ppm:
        #            tmp += i
        #            break

        tmp = 0
        for mz in notFound:
            observed.append(0)
            expected.append(isotopeSpectrum[mz])
            tmp += expected[-1]

        #print(tmp)
        expected =  expected/np.sum(expected)
        if np.sum(observed) > 0: observed = observed/np.sum(observed)

        return 1 - np.sum(np.abs(np.subtract(expected,observed))) / 2.0
        #return tmp
    except:
        return -1


def getSubForms(formula,polarity):
    f = molmass.Formula(formula)  # create formula object
    comp = f.composition()
    masses = []
    bounds = []
    carbon_pos = -1
    i = 0
    for row in comp:
        tmp = molmass.Formula(row[0])
        massC = tmp.isotope.mass
        masses.append(massC)
        if row[0] == "C":
            carbon_pos = i
        if row[0] == "H":
            if polarity == "Positive":
                bounds.append(int(row[1] + 1))
            else:
                bounds.append(int(row[1] - 1))
        else:
            bounds.append(int(row[1]))
        i += 1

    masses = np.array(masses)
    bounds = np.array(bounds)

    subs_lists = [list(range(int(b + 1))) for b in bounds]
    sub_forms = list(itertools.product(*subs_lists))

    return masses,sub_forms,carbon_pos,bounds

def createM1FromM0andFragAnnotation(spectra,masses,sub_forms,carbon_pos,bounds,ppmError):

    spectra = {str(x): val for x, val in spectra.items()}
    frags = [x for x in spectra.keys()]
    intens = [spectra[x] for x in frags]

    subsForFrags = {f: [] for f in frags}
    maxMass = -1
    for s in sub_forms:
        mass = np.dot(s, masses.transpose())
        if mass > maxMass:
            maxMass = mass
        for f in frags:
            if 1e6 * np.abs(mass - float(f)) / float(f) < ppmError:
                subsForFrags[f].append(list(s))
    subsForFrags = {key: np.array(val) for key, val in subsForFrags.items()}

    mPlus1 = {}
    unique = True
    for f, i in zip(frags, range(len(frags))):
        if len(subsForFrags[f]) > 0:
            prob_c13_in_frag = np.mean(subsForFrags[f][:, carbon_pos]) / bounds[carbon_pos]
            if len(subsForFrags[f]) > 1:
                unique = False
            mPlus1[float(f)] = intens[i] * (1 - prob_c13_in_frag)
            mz = float(f) + 1.003
            mPlus1[mz] = intens[i] * prob_c13_in_frag
    if len(mPlus1) > 0:
        maxVal = np.max(list(mPlus1.values()))
        mPlus1 = {key: val / maxVal for key, val in mPlus1.items()}

    return mPlus1,unique


def createM1SpectrumfromM0(spectra,formula,polarity,ppmError = 5):

    if len(spectra) > 0:
        try:
            masses, sub_forms, carbon_pos, bounds = getSubForms(formula,polarity)
            mPlus1,unique = createM1FromM0andFragAnnotation(spectra,masses,sub_forms,carbon_pos,bounds,ppmError)
        except:
            mPlus1 = {}
            unique = False
    else:
        mPlus1 = {}
        unique = False

    return mPlus1,unique#,{key:val[:,carbon_pos] for key,val in subsForFrags.items() if len(val) > 0}

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
                prelimClusters[x2].append([mzs[x],ids[x],groups[x]])
                good = False
                break
        if good:
            prelimClusters[mzs[x]] = [[mzs[x],ids[x],groups[x]]]
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


def getMatricesForGroup(trees,spectra,res,clusters):

    spectrum = dict(spectra[0])
    for spec in spectra[1:]:
        for m,i in spec.items():
            if m in spectrum:
                spectrum[m] += i
            else:
                spectrum[m] = i

    compoundDict = pullMostSimilarSpectra(trees,spectrum)

    keys = list(compoundDict.keys())

    if type(clusters) != type(None):
        clusterKeys = list(clusters.keys())

    metIDs = [x[0] for x in keys]
    masses = [x[1] for x in keys]
    metNames = [x[2] for x in keys]
    spectraIDs = [compoundDict[x][0] for x in keys]
    rts = [x[3] for x in keys]
    formulas = [x[4] for x in keys]
    matrix = [compoundDict[x][1] for x in keys]

    if type(clusters) != type(None):
        metIDs += [str(key[0]) + "|" + str(np.mean(key[2])) for key in clusterKeys]
        masses += [key[0] for key in clusterKeys]
        spectraIDs += [-1 for _ in clusterKeys]
        metNames += [
            "Predicted Analyte: m/z = " + str(np.round(key[0], 5)) + " rt: " +str(np.round(np.mean(key[1:]),2)) for key in clusterKeys]
        matrix += [clusters[key] for key in clusterKeys]
        rts += [np.mean(key[1:]) for key in clusterKeys]
        formulas += ["-1" for _ in clusterKeys]

    indices = list(set(flatten([list(spectrum.keys())] + [list(m.keys()) for m in matrix])))
    indices.sort()

    indicesAll = [int(np.round(x/(10**(-1*res)))) for x in indices]
    matrix = [normalizeSpectra([getVal(m, x) for x in indices]) for m in matrix]
    reduceSpec = [[getVal(spec, x) for x in indices] for spec in spectra]

    return metIDs,spectraIDs,matrix,masses,metNames,rts,formulas,indicesAll, reduceSpec

def pullMostSimilarSpectra(trees,spectra):
    returnDict = {}
    for tree,ms2Scans in trees.items():
        if len(ms2Scans) > 0:
            temp = [[id,ms2,dotProductSpectra(ms2,spectra)] for id,ms2 in ms2Scans.items()]
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

def getCompMS1Spectrum(samples,massAcc):
    # get isotope pattern match scores
    mzs = []
    for s in samples:
        if "ms1" in s:
            for mz, i in s["ms1"].items():
                mzs.append(mz)

    mzs = list(set(mzs))
    mzDict = {}
    for mz in mzs:
        bounds = [mz - massAcc * mz / 1e6, mz + massAcc * mz / 1e6]
        found = False
        for mz2 in mzDict:
            if mz2 > bounds[0] and mz2 < bounds[1]:
                found = True
                mzDict[mz2].append(mz)
                break
        if not found:
            mzDict[mz] = [mz]


    compMs1 = {}
    for mz in mzDict.keys():
        xs = []
        ys = []
        for s in samples:
            rt = s["rt"]
            ints = []
            for mz2, i in s["ms1"].items():
                if mz2 in mzDict[mz]:  # 1e6 * np.abs(mz - mz2) / mz < massAcc:
                    ints.append(i)
            if len(ints) > 0:
                xs.append(rt)
                ys.append(np.sum(ints))

        compMs1[mz] = np.trapz(ys, xs)

    return compMs1

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


def readRawDataFile(filename, maxMass, resolution, useMS1, ppmWidth = 50,offset=0.65,tic_cutoff=5,frag_cutoff=0):
   """
    Read MS datafile and convert to mzml if necessary. Conversion performs vendor centroiding. MS/MS data along with MS1
    data are extracted and returned. In addition, the contamination in the MS/MS spectra from co-isolated analytes is
    computed.

   :param filename: str, path to MS datafile
   :param maxMass: float, maximum mass to consider in MS/MS spectra
   :param resolution: int, number of decimal places to consider to m/z value of MS/MS peaks
   :param useMS1: bool, read and return MS1 data contained in file
   :param ppmWidth: float, Mass accuracy of insturment in parts per million
   :param offset: float, Isolation window width/2. Necessary for non-thermo data
   :param tic_cutoff: float, Signal cutoff for MS/MS spectra. Spectra below this signal level will be ignored
   :param frag_cutoff: float, intensity cutoff for MS/MS peaks. Fragments with absolute intensity below this threshold will be removed
   :return: result-list of MS/MS spectra. Each spectrum is a dict giving the precursor mz, retention time, isolation window
    upper and lower bounds, TIC, contamination level, and scanID
    ms1Scans-dict where key is the retention time and the value is another dict of m/z:intensity pairs

   """
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
                        try: CE = temp["precursorList"]["precursor"][0]["activation"]['collision energy']
                        except: CE = -1
                        rt = temp["scanList"]["scan"][0]["scan start time"]
                        mzs = list(zip(temp["m/z array"],temp["intensity array"]))
                        tempSpecs = []

                        spectra = {np.round(x[0],resolution):0 for x in mzs}
                        for x,y in mzs:
                            if y > frag_cutoff:
                                spectra[np.round(x,resolution)] += y

                        result.append({"id":id,"spectra":spectra,"mode":acquisitionMode,"center m/z":
                                           centerMz,"lower m/z":lowerBound,"higher m/z":upperBound,"rt":rt,"signal":tic,"CE":CE})
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
                result[samp]["ms1"] = mzScan

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


def createM1Entry(m0Entry, mode,ppm):
    key = m0Entry["id"] + "_M1"
    val = dict(m0Entry)
    val["cpdID"] = val["cpdID"] + " (M+1)"
    val["m/z"] = val["m/z"] + 1.00335
    val["name"] = val["name"] + " (M+1)"
    val["id"] = val["id"] + "_M1"
    val["spec"],_ = createM1SpectrumfromM0(m0Entry["spec"], m0Entry["formula"], mode, ppm)
    #qu.put([key, val, mode])
    return key,val,mode

class DecoID():
    """
    Class for working with and performing a database assisted deconvolution of MS/MS spectra. Deconvolution is done
    using LASSO regression.

    :param libFile: string, path to database file or "none" to use mzCloud
    :param mzCloud_lib: str, only applies to mzCloud, specifies library to seach: reference or autoprocessing
        to only use reference library
    :param numCores: int, Number of parralel processes to use.
    :param resolution: int, Number of decimal places to consider for m/z values of MS/MS peaks
    :param label: str, optional label to add to the end of output files
    :param api_key: str, for use of mzCloud api, access key must be entered
    :param mplus1PPM: float, ppm tolerance for finding subformulas in M+1 spectrum prediction. Set based on database
    :param numConcurrentGroups: int, number of unique features to processes at once, if memory consumption is high, try reducing
    :param scoringFunc: str, function to score metabolite ID matches. Defaults to the normalized dot product
    """
    def __init__(self,libFile,mzCloud_lib,numCores=1,resolution = 2,label="",api_key="none",mplus1PPM = 15,numConcurrentGroups=20,scoringFunc="dot product"):
        self.numConcurrentMzGroups = numConcurrentGroups
        self.libFile = libFile
        self.useDBLock = False
        #read tsv library file
        self.mplus1PPM = mplus1PPM
        mapper = {"dot product":dotProductSpectra,"NN":NNScoring,"shared fragment":sharedFragmentScoring}
        self.scoringFunc = mapper[scoringFunc]
        #for msp library files
        if ".msp" in libFile:
            #try:
                mz = -1
                id = -1
                cpdid = -1
                polarity = -1
                exactMass = -1
                rt = -1
                formula = -1
                spectrum = {}
                f = open(libFile,"r",encoding="utf-8")
                self.library = {"Positive":{},"Negative":{}}
                first = True
                good = False
                looking = False
                for line in f:
                    line = line.rstrip()
                    if "Name:" in line:
                        if not first:
                            if exactMass != -1 or mz != -1:
                                if polarity != -1 and formula != -1:
                                    if mz == -1:
                                        if polarity == "Negative":
                                            mz = exactMass - 1.0073
                                        else:
                                            mz = exactMass + 1.0073
                                    if cpdid == -1:
                                        cpdid = name
                                    if id != -1 and len(spectrum) > 0:
                                        self.library[polarity][id] = {"id":str(id),"cpdID":str(cpdid).replace(",","_"),"rt":rt,"formula":formula.replace(",","_"),"name":replaceAllCommas(name),"mode":polarity,"spec":spectrum,"m/z":np.round(mz,4)}
                                        good = True


                        else:
                            first = False
                        name = line.split(": ")[1]
                        mz = -1
                        id = -1
                        cpdid = -1
                        polarity = -1
                        rt = -1
                        formula = -1
                        spectrum = {}
                        looking = False
                        good = False
                    if "DB#: " in line:
                        id = line.split(": ")[1]
                    if "InChiKey:" in line or "InChIKey" in line:
                        cpdid = line.split(": ")[1]
                    if "Precursor_type:" in line:
                        temp = line.split(": ")[1]
                        if "]+" in temp:
                            polarity = "Positive"
                        elif "]+" in temp:
                            polarity = "Negative"
                    if "ExactMass:" in line:
                        exactMass = float(line.split(": ")[1])
                    if "RetentionTime:" in line:
                        rt = float(line.split(": ")[1])
                    if "Formula:" in line:
                        formula = line.split(": ")[1]
                    if "PrecursorMZ: " in line:
                        mz = float(line.split(": ")[1])
                    if "Ion_mode: " in line:
                        temp = line.split(": ")[1]
                        if temp == "P":
                            polarity = "Positive"
                        elif temp == "N":
                            polarity = "Negative"
                    if looking and len(line.split()) == 2:
                        if float(line.split()[1]) > 0:
                            spectrum[float(line.split()[0])] = float(line.split()[1])
                    if "Num Peaks" in line:
                        looking = True
                print("Library loaded successfully: " + str(len(self.library["Positive"]) + len(self.library["Negative"])) + " spectra found")
                print("Predicting M+1 now...")
                #q = Queue()
                p = Pool(numCores)
                args = []

                #processes = []
                toIts = {p:list(self.library[p].keys()) for p in self.library}
                for pol in toIts:
                    for key in toIts[pol]:
                        entry = self.library[pol][key]
                        if entry["m/z"] < maxMzForPrediction:
                            args.append([entry,pol,self.mplus1PPM])
                results = p.starmap(createM1Entry,args)
                p.close()
                p.join()
                for k,v,m in results:
                    if len(v["spec"]) > 0:
                        self.library[m][k] = v
                        #p = Process(target=createM1Entry,args=(entry,pol,q,self.mplus1PPM))
                        #while len(processes) >= numCores:
                        #    if not q.empty():
                        #        k,v,m = q.get()
                        #        if len(v) > 0:
                        #            self.library[m][k] = v
                        #    processes = [x for x in processes if x.is_alive()]
                        #p.start()
                        #processes.append(p)
                #while len(processes) > 0:
                #    if not q.empty():
                #        k, v, m = q.get()
                #        if len(v) > 0:
                #            self.library[m][k] = v
                #    processes = [x for x in processes if x.is_alive()]
                #while not q.empty():
                #    k, v, m = q.get()
                #    if len(v) > 0:
                #        self.library[m][k] = v


                pkl.dump(self.library,open(libFile.replace(".msp",".db"),"wb"))

                self.lib = customDBpy
                self.cachedReq = "none"
                self.ms_ms_library = "custom"
            # except:
            #    print(sys.exc_info())
            #    print("bad library file: ", libFile)
            #    return -1

        #for binrary datbabase files
        elif ".db" in libFile:
            self.lib = customDBpy
            self.cachedReq = "none"
            self.ms_ms_library = "custom"
            self.library = pkl.load(open(libFile,"rb"))
            print("Library loaded successfully: " + str(
                len(self.library["Positive"]) + len(self.library["Negative"])) + " spectra found")

        else:
            self.useDBLock = True
            self.lib = mzCloudPy
            self.ms_ms_library = "mzCloud"
            self.library = mzCloud_lib

        self.key = api_key
        self.numCores = numCores
        self.recursive = False
        self.label = label

    @classmethod
    def writeObj(cls,class_instance,filename):
        dill.dump(class_instance, open(filename, "wb"))

    @classmethod
    def fromDill(cls,filename):
        return dill.load(open(filename,"rb"))


    @classmethod
    def from_DecoID(cls, class_instance):
        name = str(uuid.uuid1())
        DecoID.writeObj(class_instance,name)
        obj = DecoID.fromDill(name)
        os.remove(name)
        return obj


    def readData(self,filename,resolution,peaks,DDA,massAcc,offset=.65,peakDefinitions = "",tic_cutoff=0,frag_cutoff = 0):
        """
        Read in raw MS data into DecoID object.

        :param filename: str, path to MS datafile
        :param resolution: int, number of decimal places to consider in MS/MS peaks
        :param peaks: bool, use MS1 data if available
        :param DDA: bool, True for DDA data, False for DIA data
        :param massAcc: float, mass accuracy of instrument in ppm
        :param offset: float, isolation window width measured from center of window to outer edge.
        :param peakDefinitions: str, path to peak definition file
        :param tic_cutoff: float, TIC cutoff for MS/MS spectra. Spectra below this value will be ignored
        :param frag_cutoff: float, Intensity cutoff for fragments to be included. Fragments in spectra below this intensity will be ignored.
        :return: None
        """
        #read in data
        samples, ms1 = readRawDataFile(filename, MAXMASS, resolution, peaks,offset=offset,tic_cutoff=tic_cutoff,ppmWidth=massAcc,frag_cutoff=frag_cutoff)
        if samples != -1:
            #set object fields
            self.samples = samples
            self.ms1 = ms1
            self.resolution = resolution
            self.peaks = peaks
            self.massAcc = massAcc
            self.DDA = DDA
            if len(self.ms1) < 1:
                self.peaks = False
            if not DDA and not peaks:
                self.DDA = True
            # samples = samples[:100]
            #output # of spectra collected
            print(len(self.samples), " MS2 spectra detected")
            fileending = "." + filename.split(".")[-1]
            self.filename = filename.replace(fileending, "")
            self.filename = self.filename.replace('"', "")
            self.clusters = {}
            self.data2WriteFinal = []
            self.outputDataFinal = {}

            #if peak definitions are provided
            if peakDefinitions != "":
                self.peakData = pd.read_csv(peakDefinitions)
                self.peakData = self.peakData[["mz","rt_start","rt_end"]] #read peak info
                dfIndex = list(self.peakData.index.values)

                if self.DDA: #if DDA group spectra using peak information
                    i = -1
                    goodSamps = []
                    self.samples = []
                    for index in dfIndex:
                        for samp in samples:
                            if abs(samp["center m/z"]-self.peakData.at[index,"mz"])/self.peakData.at[index,"mz"]*1e6 < self.massAcc:
                                if samp["rt"] >= self.peakData.at[index,"rt_start"] and samp["rt"] <= self.peakData.at[index,"rt_end"]:
                                    newSamp = dict(samp)
                                    newSamp["group"] = index
                                    self.samples.append(newSamp)

                    numGroups = len(set([samp["group"] for samp in self.samples]))
                    print("Number of compounds with acquired MS2: ",numGroups)
                    print("Number of spectra to deconvolve: ",len(self.samples))

                #if DIA
                else:
                    samplesDict = {x["id"]: x for x in samples} #get unique windows
                    mzGroups = clusterSpectraByMzs([x["center m/z"] for x in samples], [x["id"] for x in samples],
                                                   [x["rt"] for x in samples],list(range(len(samples))) ,self.massAcc)
                    index = 0 #put samples into groups
                    for mz in mzGroups:
                        for m, id,g in mzGroups[mz]:
                            samplesDict[id]["group"] = index
                        index += 1
            #if peak information is not provided set each ms/ms spectrum as a unique feature
            else:
                index = 0
                for samp in samples:
                    samp["group"] = index
                    index += 1
                self.peakData = pd.DataFrame.from_dict({i:{"mz":samp["center m/z"],"rt_start":samp["rt"]-.05,"rt_end":samp["rt"]+.05} for i,samp in zip(range(len(samples)),samples)},orient="index")


    def readMS_DIAL_data(self,file,mode,massAcc,peakDataFile):
        """
        Load in data from MS-DIAL exported peak list text file of deconvoluted spectra. This enables library searching and
        combined usage of DecoID and MS-DIAL

        :param file: str, path to MS-DIAL output file
        :param mode: str, polarity (Positive/Negative)
        :param massAcc: int, mass accuracy of instrument (ppm)
        :param peakDataFile: str, path to peak information file for features of interest. If doing a combined usage, this needs to be the same file used for DecoID>
        :return:
        """
        self.resolution = 2
        self.peaks = False
        DDA = True
        self.peakData = pd.read_csv(
           peakDataFile)
        self.peakData = self.peakData[["mz", "rt_start", "rt_end"]]
        i = -1
        sampleData = pd.read_csv(file, sep="\t")
        samples = []

        for index, row in sampleData.iterrows():
            id = row["PeakID"]
            if not pd.isna(row["MSMS spectrum"]):
                spectra = row["MSMS spectrum"].split()
                spectra = [x.split(":") for x in spectra]
                spectra = {float(mz): float(i) for mz, i in spectra if float(i) > 0}
                spectra = collapseAsNeeded(spectra, 2)
                centerMz = row["Precursor m/z"]
                upperBound = centerMz + .1
                lowerBound = centerMz - .1
                rt = row["RT (min)"]
                signal = row["Area"]
                samples.append({"id": id, "spectra": spectra, "mode": mode, "center m/z":
                    centerMz, "lower m/z": lowerBound, "higher m/z": upperBound, "rt": rt, "signal": signal})

        self.ms1 = {}
        self.massAcc = massAcc
        self.DDA = DDA

        # samples = samples[:100]
        fileending = "." + file.split(".")[-1]
        self.filename = file.replace(fileending, "")
        self.filename = self.filename.replace('"', "")
        self.clusters = {}
        self.data2WriteFinal = []
        self.outputDataFinal = {}

        goodSamps = []
        for samp in samples:
            i += 1
            for index, row in self.peakData.iterrows():
                if abs(samp["center m/z"] - row["mz"]) / row["mz"] * 1e6 < massAcc and samp["rt"] >= row[
                    "rt_start"] and samp["rt"] <= row["rt_end"]:
                    samp["group"] = index
                    goodSamps.append(i)
                    break
        self.samples = [samples[x] for x in goodSamps]
        print(len(self.samples), " MS2 spectra detected")

    def readDecoMS2_data(self, file, mode, massAcc, peakDataFile):
        self.resolution = 2
        self.peaks = False
        DDA = True
        self.peakData = pd.read_csv(
            peakDataFile)
        self.peakData = self.peakData[["mz", "rt_start", "rt_end"]]
        i = -1
        sampleData = pd.read_csv(file)
        samples = []

        for index, row in sampleData.iterrows():
            id = row["index"]
            spectra = row["spectrum"].split()
            spectra = [x.split(":") for x in spectra]
            spectra = {float(mz): float(i) for mz, i in spectra if float(i) > 0}
            spectra = collapseAsNeeded(spectra, self.resolution)
            centerMz = row["mz"]
            upperBound = centerMz + .1
            lowerBound = centerMz - .1
            rt = (row["rt_start"]+row["rt_end"])/2
            signal = 1e6
            samples.append({"id": id, "spectra": spectra, "mode": mode, "center m/z":
                centerMz, "lower m/z": lowerBound, "higher m/z": upperBound, "rt": rt, "signal": signal})

        self.ms1 = {}
        self.massAcc = massAcc
        self.DDA = DDA

        # samples = samples[:100]
        fileending = "." + file.split(".")[-1]
        self.filename = file.replace(fileending, "")
        self.filename = self.filename.replace('"', "")
        self.clusters = {}
        self.data2WriteFinal = []
        self.outputDataFinal = {}

        goodSamps = []
        for samp in samples:
            i += 1
            for index, row in self.peakData.iterrows():
                if abs(samp["center m/z"] - row["mz"]) / row["mz"] * 1e6 < massAcc and samp["rt"] >= row[
                    "rt_start"] and samp["rt"] <= row["rt_end"]:
                    samp["group"] = index
                    goodSamps.append(i)
                    break
        self.samples = [samples[x] for x in goodSamps]
        print(len(self.samples), " MS2 spectra detected")

    def recursiveRunQ(self,q, data2WriteLater,outputDataFile):
        toQuit = False
        while True:
            if not q.empty():
                new = q.get()

                if new == "EXIT" or toQuit:
                    toQuit = True
                    if q.empty():
                        return 0
                [result, centerMz, id, rt, s2n, indices, numComponents, fragments, decoSpec,components] = new
                if fragments == "none":
                    fragments = [centerMz]
                data2WriteLater.append(list(new))
                outputDataFile[id] = {}
                outputDataFile[id]["decoSpec"] = decoSpec
                outputDataFile[id]["center m/z"] = centerMz
                outputDataFile[id]["rt"] = rt
                outputDataFile[id]["Hits"] = {}
                outputDataFile[id]["index"] = indices
                outputDataFile[id]["fragments"] = fragments
                for x in result:
                    outputDataFile[id]["Hits"][x] = result[x]
            else:
                time.sleep(1)

    def identifyUnknowns(self,resPenalty=100,percentPeaks=0.01,iso=False,ppmThresh = 10,dpThresh = 80,rtTol=.5):
        """
        Generate the on-the-fly unknown library by searching all spectra first and identifying those that are unknown.
        These spectra can then be used to deconvolve other spectra. Only applicable to DDA with MS1 collected.

        :param resPenalty: float, lasso regression coefficient. Set to float("inf") for direct searching. Set to 0 for unregularized deconvolution. Reccomend 1 for DDA, 100 for DIA
        :param percentPeaks: float, filtering parameter, exclude spectra that do not match this fraction of the database peaks
        :param iso: bool, remove contamination from orphan isotopologues if True, False- do not.
        :param ppmThresh: float, mass error match tolerance for the resulting hits of the spectra.
        :param dpThresh: float, dot product match tolerance for the resulting hits of the spectra.
        :param rtTol: float, retention time tolerance (minutes) allowed between database and acquired RT in order for database spectrum to be used in deconvolution
        :return: None
        """
        try: self.samples
        except: print("datafile not loaded");return -1
        self.iso = iso
        self.percentPeaks = percentPeaks
        self.rtTol = rtTol
        self.resPenalty = resPenalty
        self.redundancyCheckThresh = np.inf

        #check datatype
        q = Queue(maxsize=1000)
        if not self.peaks or not self.DDA:
            print("MS1 and DDA Required for Unknown ID Recursion")
            self.recursive = False
            return -1
        else:
            self.recursive = True
            samplesToGo = []
            samplesLeft = []
            sigThresh = 4
            for x in self.samples:
                if (x["percentContamination"] < .1) and x["signal"] > sigThresh: #identify non-chimeric spectra with a TIC > 1e4
                    samplesToGo.append(x)
                else:
                    samplesLeft.append(x)

            data2Write = []
            outputData = {}
            t = Thread(target=self.recursiveRunQ, args=(q,data2Write,outputData))
            t.start()
            #search spectra in questions
            self.runSamples(samplesToGo,q)

            t.join()

            #filter results for known spectra
            unknownGroupIDs = []
            unknownThreshold = dpThresh
            unknownPPMThreshold = ppmThresh
            for result, centerMz, id, rt, s2n, indices, numComponents, fragments, decoSpec,components in data2Write:
                if numComponents < 4:
                    if len(result) > 0:
                        if not any(result[x][0] >= unknownThreshold and (abs(float(result[x][4]) - float(x[2])) * 1e6) / float(x[2]) < unknownPPMThreshold for x in result):
                            unknownGroupIDs.append(id)
                    else:
                        unknownGroupIDs.append(id)

            unknownSamples = []
            for id in unknownGroupIDs:
                for samp in samplesToGo:
                    if samp["group"] == id:
                        unknownSamples.append(samp)

            #add spectra to unknown library for future usage
            self.clusters = clusterSpectra(unknownSamples,self.peakData)

            print("Number of Predicted Unknown Compounds: ", len(self.clusters))



    def searchSpectra(self,verbose,resPenalty = 100,percentPeaks=0.01,iso=False,threshold = 0.0,rtTol=.5,redundancyCheckThresh = 0.9):
        """
        Search the spectra loaded into the DecoID object and write the output files

        :param verbose: str or Queue, "y" to write the progress to std out. Queue is used with the GUI to send updates dynamically/
        :param resPenalty: float, lasso regression coefficient. Set to float("inf") for direct searching. Set to 0 for unregularized deconvolution. Reccomend 1 for DDA, 100 for DIA
        :param percentPeaks: float, filtering parameter, exclude spectra that do not match this fraction of the database peaks
        :param iso: bool, remove contamination from orphan isotopologues if True, False- do not.
        :param threshold: float filtering parameter to remove hits with a spectral similarity less than this value
        :param rtTol: float, retention time tolerance to use database spectra (minutes)
        :param redundancyCheckThresh: float, dot product threshold to classify a component as redundant set to > 100 to turn off
        :return: None
        """
        try: self.samples
        except: print("datafile not loaded");return -1
        q = Queue(maxsize=1000)

        self.iso = iso
        self.percentPeaks = percentPeaks
        self.rtTol = rtTol
        self.redundancyCheckThresh = redundancyCheckThresh

        self.resPenalty = resPenalty
        if len(self.ms1) < 1:
            self.recursive = False
            self.peaks = False
            # if resPenalty > 10:
            #     self.resPenalty = 10

        if not self.peaks and not self.DDA:
            self.DDA = True


        if type(self.samples) != type(-1):

            success = False
            error = True
            #check output file can be opened.
            while (not success):
                try:
                    outfile = open(self.filename+ self.label + "_decoID" + ".csv", "w",encoding="utf-8")
                    success = True
                except:
                    if error:
                        print("Error: please close " + self.filename + self.label + "_decoID.csv")
                        error = False
                        time.sleep(2)
            #write header
            outfile.write(
                "#featureID,isolation_center_m/z,rt,compound_m/z,compound_rt,compound_formula,DB_Compound_ID,Compound_Name,DB_Spectrum_ID,dot_product,ppm_Error,isotope_similiarity,Abundance,ComponentID,redundant\n")
            #start Q
            t = Thread(target=self.runQ, args=(
            q,outfile,threshold,self.outputDataFinal,self.filename + self.label + ".DecoID",
            self.filename + self.label + "_scanInfo.csv", verbose))
            t.start()

            t2 = Thread(target=self.sendCompleteSamples, args=(self.data2WriteFinal,q))
            t2.start()
            #search spectra
            self.runSamples(self.samples,q)

            #wait for completion
            t.join()
            t2.join()


    def sendCompleteSamples(self, data2Send,q):
        for x in data2Send:
            success = False
            while not success:
                try:
                    q.put(x, timeout=5)
                    success = True
                except:
                    pass

    def runQ(self,q,outfile,threshold, outputDataFile, filename, scanInfoFileName,
             verbose="y"):
        index = 0
        status = 0
        outputScanFile = open(scanInfoFileName, "w",encoding="utf-8")
        outputScanFile.write("#featureID,Signal to Noise Ratio,numComponents,componentID,componentRT,componentAbundance,componentMz,spectrum\n")
        toQuit = False
        delimiter = ","
        if self.DDA:
            numSamples = len(set([samp["group"] for samp in self.samples]))
        else:
            numSamples = len(self.peakData)

        indexRef = [np.round(x * 10 ** (-1 * self.resolution), self.resolution) for x in list(range(int(MAXMASS * 10 ** self.resolution)))]

        while True:
            if not q.empty():
                new = q.get()
                if new == "EXIT" or toQuit:
                    toQuit = True
                    if q.empty():
                        outfile.close()
                        if filename != "None":
                            fh = gzip.open(filename, "wb")
                            pkl.dump([self.samples,self.peakData,self.ms1,outputDataFile], fh)
                            fh.close()
                        outputScanFile.close()
                        break
                elif new != "EXIT":
                    [result, centerMz, id, rt, s2n, indices, numComponents, fragments, decoSpec,components] = new

                    if fragments == "none":
                        fragments = [centerMz]

                    for c in components:
                        outputScanFile.write(str(id) + "," + str(s2n) + "," + str(numComponents) + "," + str(c[0]) + "," +str(c[2]) + "," + str(c[3]) + "," + str(c[1]) + ",")
                        for ind,i in zip(indices,components[c]):
                            if i > 0:
                                outputScanFile.write(str(np.round(ind * 10 ** (-1 * self.resolution), self.resolution)) + ":" + str(i) + " ")
                        outputScanFile.write("\n")
                    update = False
                    if id not in outputDataFile:
                        outputDataFile[id] = {}
                        outputDataFile[id]["NumCand"] = len(result)
                        outputDataFile[id]["decoSpec"] = decoSpec
                        outputDataFile[id]["center m/z"] = centerMz
                        outputDataFile[id]["rt"] = rt
                        outputDataFile[id]["Hits"] = {}
                        outputDataFile[id]["index"] = indices
                        outputDataFile[id]["fragments"] = fragments
                        outputDataFile[id]["components"] = components
                        update = True
                    for x in result:
                        if result[x][0] > threshold:
                            try:
                                tempMz = result[x][4]
                            except:
                                print(result[x])
                                break

                            componentName = result[x][3]
                            outfile.write(str(id) + "," + str(tempMz))
                            [outfile.write(delimiter + z) for z in
                             [str(rt), str(x[2]),str(x[4]),str(x[5]), str(x[0]), str(x[3]), str(x[1]), str(result[x][0]),
                              str((float(x[2]) - float(tempMz)) * 1e6 / float(tempMz)), str(x[6]),str(x[7]),str(componentName),str(result[x][6])]]

                            outfile.write("\n")
                            if filename != "None" and update:
                                outputDataFile[id]["Hits"][x] = result[x]
                index += 1
                if type(verbose) == type(str()):
                    if "y" in verbose:
                        if status == 0:
                            print("0................................................100")
                        if numSamples == 0:
                            for _ in range(50):
                                if "yV" in verbose:
                                    print("x", flush=True)
                                else:
                                    print("x", end="", flush=True)
                        else:
                            while index / numSamples > status:
                                if "yV" in verbose:
                                    print("x", flush=True)
                                else:
                                    print("x", end="", flush=True)
                                status += .02

                else:
                    if numSamples == 0:
                        for _ in range(50):
                            verbose.put("1")
                    else:
                        while index / numSamples > status:
                            status += 0.02
                            verbose.put("1")

            else:
                time.sleep(1)

    def prepareForCluster(self,numberOfFiles):
        """
        Split a large dataset read into DecoID for searching on separate machines. Helpful for large datafiles that want to be run on a compute cluster.
        Will create dill files of DecoID objects with partial data lists, but the same parameters. These can be read with fromDill

        :param numberOfFiles: int, number of data files to split
        :return: None
        """
        # if self.recursive:
        #     tempObject = DecoID.from_DecoID(self)
        #     tempObject.filename+="_unknowns_"
        #     tempObject.samples = []
        #     DecoID.writeObj(tempObject,tempObject.filename+".dill")
        #get number of groups
        groups = list(set([x["group"] for x in self.samples]))
        numGroupsPerFile = int(len(groups)/numberOfFiles)
        self.data2WriteFinal = []
        self.outputDataFinal = {}
        groupBegin = 0
        processes = []
        groupEnd = groupBegin + numGroupsPerFile
        for num in range(numberOfFiles): #for each file
            tempObject = DecoID.from_DecoID(self)
            tempObject.filename += "_" + str(num) + "_"
            if num != numberOfFiles - 1: #specify groups
                group = groups[groupBegin:groupEnd]
            else:
                group = groups[groupBegin:]
            groupBegin = groupEnd
            groupEnd = groupBegin + numGroupsPerFile #fix samples
            tempObject.samples = [samp for samp in tempObject.samples if samp["group"] in group]
            DecoID.writeObj(tempObject,tempObject.filename+".dill") #write file


    @staticmethod
    def combineResultsAutoAnnotate(filenames, newFilename, numHits=1,min_score=0):
        """
        Combine results from several files. Takes the best hit found across multiple files

        :param filenames: list, filenames to merge should not have any file endings (i.e. .mzML, _decoID.csv)
        :param newFilename: output file name
        :param numHits: int, number of hits to return in the merged file for each feature
        :param min_score: float, minimum dot product score for returned hits.
        :return: None
        """
        bestResults = pd.read_csv(filenames[0] + "_decoID.csv")
        bestResults = bestResults[bestResults["DB_Spectrum_ID"] != "-1"]
        bestResults = bestResults[bestResults["DB_Spectrum_ID"] != -1]
        bestResults = bestResults[bestResults["dot_product"] > min_score]
        bestResults["Origin File"] = [filenames[0] for x in range(len(bestResults))]

        feats = list(set(bestResults["#featureID"].values))
        for f in filenames[1:]:
            data = pd.read_csv(f + "_decoID.csv")
            data = data[data["DB_Spectrum_ID"] != "-1"]
            data = data[data["DB_Spectrum_ID"] != -1]
            data = data[data["dot_product"] > min_score]
            data["Origin File"] = [f for x in range(len(data))]

            found_feats = list(set(data["#featureID"].values))
            for feat in found_feats:
                part = data[data["#featureID"] == feat]
                if feat not in feats:
                    bestResults = pd.concat((bestResults, part), axis=0, ignore_index=True)
                    feats.append(feat)
                else:
                    tmp = bestResults[bestResults["#featureID"] == feat]
                    toAdd = []
                    for index, row in part.iterrows():
                        if row["DB_Compound_ID"] in tmp["DB_Compound_ID"].values:
                            score = [[i, r["dot_product"]] for i, r in tmp.iterrows() if
                                     r["DB_Compound_ID"] == row["DB_Compound_ID"]][0]
                            if row["dot_product"] > score[1]:
                                bestResults.loc[score[0], :] = row
                        else:
                            toAdd.append(index)
                    if len(toAdd) > 0:
                        bestResults = pd.concat((bestResults,part.loc[toAdd,:]),axis=0,ignore_index=True)

        toDrop = []
        for feat in feats:
            tmp = bestResults[bestResults["#featureID"] == feat]
            tmp = tmp.sort_values(by="dot_product", ascending=False)
            if len(tmp) > numHits:
                toDrop += list(tmp.index.values[numHits:])

        bestResults = bestResults.drop(toDrop)
        bestResults = bestResults.sort_values(by="#featureID", axis=0)
        bestResults.to_csv(newFilename + "_decoID.csv", index=False)

    @staticmethod
    def combineResults(filenames,newFilename,endings = ["_scanInfo.csv","_decoID.csv",".DecoID"]):
        """
        Combine results from several files. Helpful to merge results from lots of files. Does only the .csv files for memory reasons.

        :param filenames: filenames to merge
        :param newFilename: output file name
        :param endings: which file endings to use.
        :return: None
        """
        goodFiles = []
        for file in filenames:
            try:
                open(file+"_decoID.csv","r")
                goodFiles.append(file)
            except:
                pass
        if len(goodFiles) > 0:
            for end in endings:
                if ".csv" in end:
                    base = open(goodFiles[0]+end,"r").readlines()
                    if len(goodFiles) > 1:
                        for f in goodFiles[1:]:
                            base += open(f+end,"r").readlines()[1:]
                        outfile = open(newFilename+end,"w")
                        [outfile.write(x) for x in base]
                        outfile.close()
        else:
            print("no files detected")
    def runSamples(self,samples,q):
        numConcurrentMzGroups = self.numConcurrentMzGroups
        numProcesses = Value('i', 0)
        lock = Lock()
        processes = []
        dbLock = Lock()
        availableToGrab = Value('i',1)
        samplesDict = {(x["id"],x["group"]):i for x,i in zip(samples,range(len(samples)))}
        keys = list(samplesDict.keys())
        mzGroups = clusterSpectraByMzs([samples[samplesDict[k]]["center m/z"] for k in keys],[samples[samplesDict[k]]["id"] for k in keys],[samples[samplesDict[k]]["rt"] for k in keys],[samples[samplesDict[k]]["group"] for k in keys],self.massAcc)
        mzGroupsWSamples = {x:{"samples":[]} for x in mzGroups}
        for mz in mzGroups:
            for m,id,g in mzGroups[mz]:
                mzGroupsWSamples[mz]["samples"].append(samples[samplesDict[(id,g)]])

        for mz in mzGroupsWSamples:
            mzGroupsWSamples[mz]["lower m/z"] = np.min([x["lower m/z"] for x in mzGroupsWSamples[mz]["samples"]])
            mzGroupsWSamples[mz]["higher m/z"] = np.max([x["higher m/z"] for x in mzGroupsWSamples[mz]["samples"]])
            if type(self.clusters) != type(None):
                toAdd = [c for c in self.clusters if
                         c["m/z"] >= mzGroupsWSamples[mz]["lower m/z"] and c["m/z"] <= mzGroupsWSamples[mz]["higher m/z"]]
                toAdd = {(c["m/z"], c["rtWindow"][0], c["rtWindow"][1]): c["spectrum"] for
                         c in toAdd}
            else:
                toAdd = self.clusters

            #self.processMzGroup(mzGroupsWSamples[mz]["samples"], toAdd, numProcesses, mzGroupsWSamples[mz]["higher m/z"],
            #    mzGroupsWSamples[mz]["lower m/z"], availableToGrab, lock, q)
            p = Thread(target=self.processMzGroup,args=(mzGroupsWSamples[mz]["samples"],toAdd,numProcesses,mzGroupsWSamples[mz]["higher m/z"],mzGroupsWSamples[mz]["lower m/z"],availableToGrab,lock,q,dbLock))
            while len(processes) >= numConcurrentMzGroups:
                processes = [proc for proc in processes if proc.is_alive()]
                time.sleep(1)
            p.start()
            processes.append(p)
        for p in processes:
            p.join()

        q.put("EXIT")


    def processMzGroup(self,allSamples, clusters, numProcesses, upperBound, lowerBound, availableToGrab,lock,q,dbLock):

        def startProc(numP, l,samples,trees):
            spectra = [samp["spectra"] for samp in samples]
            toAdd = {key: value for key, value in clusters.items() if
                     any(key[1] < samp["rt"] and key[2] > samp["rt"] for samp in samples)}
            p = Process(target=DecoID.processGroup,
                        args=(
                            samples,group,trees, spectra, self.iso,self.DDA,self.massAcc,self.peaks,q,self.resPenalty,self.resolution,toAdd,self.peakData,lowerBound,upperBound,self.rtTol,self.redundancyCheckThresh,self.scoringFunc))
            #t = Thread(target=startProc, args=(p, numProcesses, lock))
            while not checkRoom(numP, l):
                time.sleep(1)
            p.start()
            p.join()
            with l:
                numP.value -= 1
            return 0


        def checkRoom(val, l):
            with l:
                if val.value < self.numCores:
                    #print(val.value,self.numCores)
                    val.value += 1
                    return True
                else:
                    return False

        mode = allSamples[0]["mode"]
        if self.useDBLock:
            with dbLock:
                trees = self.lib.getCanidateCompoundTrees(mode, upperBound, lowerBound, self.iso,
                                                                          self.library,self.key,self.massAcc,self.resolution)
        else:
            trees = self.lib.getCanidateCompoundTrees(mode, upperBound, lowerBound,
                                                                                   self.iso,
                                                                                   self.library, self.key,self.massAcc,self.resolution)
        featureGroups = {x["group"]:[] for x in allSamples}
        [featureGroups[x["group"]].append(x) for x in allSamples]



        threads = []
        for group,samples in featureGroups.items():
            t = Thread(target=startProc,args=(numProcesses,lock,samples,trees))
            t.start()
            threads.append(t)

        for t in threads:
            t.join()



    @staticmethod
    def processGroup(samples,group,trees, spectra, iso,DDA,massAcc,peaks,q,resPenalty,resolution,toAdd,peakData,lowerbound,upperbound,rtTol,redundancyCheckThresh,scoringFunc):
        [metIDs, spectraIDs, matrix, masses, metNames, rts, formulas, indicesAll, reduceSpec] = getMatricesForGroup(trees, spectra,resolution, toAdd)
        isoIndices = [x for x in range(len(metIDs)) if "(M+1)" in metNames[x]]



        results = []
        samples2Go = list(samples)
        reduceSpec2Go = list(reduceSpec)

        for sample, spectrum in zip(samples2Go, reduceSpec2Go):
            results.append(DecoID.processSample(sample, spectrum, masses, matrix, massAcc, peaks,resPenalty,isoIndices,rts,rtTol))

        if DDA:
            combinedSpectrum = np.sum(reduceSpec, axis=0)
            resVector = np.sum([x[0] for x in results], axis=0)
            centerMz = np.mean([s["center m/z"] for s in samples])
            rt = np.mean([s["rt"] for s in samples])
            foundSpectra = np.sum([x[2] for x in results], axis=0)
            numComponents = len([i for i in resVector if i > 1e-8])
            s2n = (1.0 / max([1, numComponents])) * np.sum(flatten(foundSpectra)) / np.sum(
                np.abs(np.subtract(foundSpectra, combinedSpectrum)))
            if peaks:
                frags = flatten([s["fragments"] for s in samples])
            else:
                frags = []

            scores, components = scoreDeconvolution(combinedSpectrum, matrix, resVector, metIDs,
                                                    masses, centerMz,rts, massAcc,rtTol,rt,redundancyCheckThresh,scoringFunc,indicesAll,resolution,sample["mode"])

            compMs1 = getCompMS1Spectrum(samples,massAcc)

            isotopeScores = [scoreIsotopePattern(compMs1, f, samples[0]["mode"], massAcc) for f in formulas]

            result = {(met, specID, mass, name,r,formula,isoScore, safeDivide(comp, sum(resVector))): coeff for
                      met, name, specID, coeff, mass, comp,r,formula,isoScore in
                      zip(metIDs, metNames, spectraIDs, scores, masses, resVector,rts,formulas,isotopeScores)}

            # output result
            success = False
            while not success:
                try:
                    q.put([result, centerMz, group, rt, s2n, indicesAll, len(components), frags, foundSpectra, components],
                          timeout=5)
                    success = True
                except:
                    print("waiting to put in q")
                    pass
        else:
            rts1 = [samp["rt"] for samp in samples]
            for index,row in peakData.iterrows():
                if row["mz"] >= lowerbound and row["mz"] <= upperbound:
                    goodIndices = [x for x in range(len(rts1)) if rts1[x] >= row["rt_start"] and rts1[x] <= row["rt_end"]]
                    if len(goodIndices) > 0:
                        combinedSpectrum = np.sum([reduceSpec[x] for x in goodIndices], axis=0)
                        resVector = np.sum([results[x][0] for x in goodIndices], axis=0)
                        centerMz = row["mz"]
                        rt = np.mean([row["rt_start"], row["rt_end"]])
                        foundSpectra = np.sum([results[x][2] for x in goodIndices], axis=0)
                        numComponents = len([i for i in resVector if i > 1e-8])
                        s2n = (1.0 / max([1, numComponents])) * np.sum(flatten(foundSpectra)) / np.sum(
                            np.abs(np.subtract(foundSpectra, combinedSpectrum)))
                        if peaks:
                            frags = flatten([s["fragments"] for s in samples])
                        else:
                            frags = []

                        scores, components = scoreDeconvolution(combinedSpectrum, matrix, resVector, metIDs,
                                                                masses, centerMz, rts, massAcc, rtTol, rt,redundancyCheckThresh,scoringFunc,indicesAll,resolution,sample["mode"])

                        compMs1 = getCompMS1Spectrum([samples[x] for x in goodIndices], massAcc)

                        isotopeScores = [scoreIsotopePattern(compMs1, f, samples[0]["mode"], massAcc) for f in formulas]

                        result = {(met, specID, mass, name, r, formula,isoScore, safeDivide(comp, sum(resVector))): coeff for
                                  met, name, specID, coeff, mass, comp, r, formula,isoScore in
                                  zip(metIDs, metNames, spectraIDs, scores, masses, resVector, rts, formulas,isotopeScores)}

                        # output result
                        success = False
                        while not success:
                            try:
                                q.put([result, centerMz, index, rt, s2n, indicesAll, len(components), frags, foundSpectra,
                                       components],
                                      timeout=5)
                                success = True
                            except:
                                print("waiting to put in q")
                                pass

    @staticmethod
    def processSample(sample,spectrum,masses,matrix,massAcc,peaks,resPenalty,isoIndices,rts,rtTol):

        centerMz = sample["center m/z"]
        lowerBound = sample["lower m/z"]  # lowest m/z in isolation window
        upperBound = sample["higher m/z"]  # highest
        rt = sample["rt"]
        def checkScan(mass,index):
            if index in isoIndices:
                return inScanIso(sample["fragments"],mass,massAcc)
            else:
                return inScan(sample["fragments"],mass,massAcc)

        if len(matrix) > 0:
            goodIndices = [x for x in range(len(masses)) if masses[x] < upperBound and masses[x] > lowerBound and ((not peaks) or checkScan(masses[x],x)) and (np.abs(rts[x]-rt) < rtTol or rts[x] == -1)]
            if len(goodIndices) > 0:
                resTemp, s2n = solveSystem([matrix[x] for x in goodIndices], spectrum,resPenalty)  # deconvolution
                res = [0.0 for _ in masses]
                for x in range(len(goodIndices)):
                    res[goodIndices[x]] = resTemp[x]
                decoSpec = flatten(np.dot(np.transpose(matrix), [[x] for x in res]).tolist())
            else:
                res = [0.0 for _ in masses]
                s2n = 0.0
                decoSpec = [0 for _ in spectrum]
            return [res,s2n,decoSpec]
        else:
            return [[],0,[0 for _ in spectrum]]

    @staticmethod
    def toSiriusOutput(filename,polarity,spectype="o",ppmErr=10):

        data = pd.read_csv(filename)
        dir = filename.replace("_scanInfo.csv","_sirius/")
        os.mkdir(dir)
        if spectype == "o":
            key = "original"
            data = data[data["componentID"] == key]
        scanIDs = list(set(data["#scanID"].values))
        for scanID in scanIDs:
            relevant = data[data["#scanID"] == scanID]
            targetMz = relevant[relevant["componentID"] == "original"]
            targetMz = targetMz.at[targetMz.index.values[0],"componentMz"]
            for index,row in relevant.iterrows():
                if 1e6*abs(row["componentMz"] - targetMz)/targetMz < ppmErr:
                    spec = row["spectrum"]
                    id  = row["componentID"].replace("|","-")
                    id = id.replace("/","-")
                    outfile = open(dir+str(scanID) + "_" + str(id) + ".ms","w")
                    outfile.write(">compound " + str(scanID) + "_" + str(id))
                    outfile.write("\n>parentmass " + str(row["componentMz"]))
                    outfile.write("\n>charge " + str(polarity))
                    outfile.write("\n>ms1\n")
                    outfile.write(str(row["componentMz"]) + " 100")
                    outfile.write("\n\n\n>ms2")
                    for pair in spec.split():
                        mz = pair.split(":")[0]
                        i = pair.split(":")[1]
                        outfile.write("\n" + mz + " " + i)
                    outfile.close()



class customDBpy():
    def __init__(self):
        pass
    @staticmethod
    def getCanidateCompoundTrees(mode, upperBound, lowerBound, isotope=False, library="none", key="none",
                                 ppmError=5,res=2):

        possIsotopes = {}
        # get compounds in isolation window
        possCompounds = {
            (library[mode][x]["cpdID"],
             library[mode][x]["m/z"],library[mode][x]["name"],
             library[mode][x]["rt"],library[mode][x]["formula"],
             library[mode][x]["id"]):library[mode][x]["spec"] for x in library[mode] if library[mode][x]["m/z"] >= lowerBound and library[mode][x]["m/z"] <= upperBound and "(M+1)" not in library[mode][x]["name"]}


        # get isotopes if necessary
        if isotope:
            possIsotopes = {(library[mode][x]["cpdID"],
                             library[mode][x]["m/z"], library[mode][x]["name"],
                             library[mode][x]["rt"],library[mode][x]["formula"],
                             library[mode][x]["id"]):library[mode][x]["spec"] for x in library[mode] if library[mode][x]["m/z"]-1.0035 >= lowerBound - 1.00335 and library[mode][x]["m/z"]-1.0035 <= lowerBound and "(M+1)" in library[mode][x]["name"]}

        possCompounds.update(possIsotopes)
        cpds = []
        trees = {}
        for x in possCompounds:
            if x[0] not in cpds:
                cpds.append(x[0])
                trees[x[:-1]] = {}
                id = x[:-1]
            else:
                id = [t for t in trees if t[0] == x[0]][0]

            trees[id][x[-1]] = collapseAsNeeded(possCompounds[x],res)

        return trees

class Keys():
    def __init__(self,api_code):
        self.HEADER = {
            'V1-API-Secret': api_code,
            'cache-control': "no-cache",
            }
    SPECTRAURL = "https://mzcloud.org/api/v1/LIBRARY/trees/TREENUMBER/spectra"
    COMPOUNDURL = "https://mzcloud.org/api/v1/LIBRARY/compounds"

class mzCloudPy():
    def __init__(self):
        pass
    @staticmethod
    def getCanidateCompoundTrees(mode, upperBound, lowerBound, isotope=False, library="reference", key="none",ppmError=5,res=2):

        keys = Keys(key)
        # make list of isolation window range
        possIsotopes = set()
        possCompounds = set()

            # get compounds in isolation window
        possCompounds = possCompounds.union(set([tuple(x) for x in MZCOMPOUNDTREELINK[library][mode] if
                             float(x[5]) >= lowerBound and float(x[5]) <= upperBound]))

        # get isotopes if necessary
        if isotope:
            possIsotopes = possIsotopes.union(set([tuple(x) for x in MZCOMPOUNDTREELINK[library][mode] if
             float(x[5]) >= lowerBound - 1.00335 and float(x[5]) <= lowerBound]))

        possCompounds = possCompounds.union(possIsotopes)
        trees = mzCloudPy.getTrees(possCompounds, keys,library=library)
        possIsotopes = list(possIsotopes)
        returnDict = {}
        cpds = []
        for tree in trees:
            if tree in possIsotopes:
                name = tree[2] + " (M+1)"
                mz = tree[5] + 1.00335
                id = library[0] + str(tree[1]) + " (M+1)"
            else:
                name = tree[2]
                mz = tree[5]
                id = library[0] + str(tree[1])
            if id not in cpds:
                cpds.append(id)
                returnDict[(id,mz,name,tree[4],tree[3])] = {}

            if tree in possIsotopes:
                masses, sub_forms, carbon_pos, bounds = getSubForms(tree[3],mode)
            for specID in trees[tree]:
                if tree in possIsotopes:
                    spec = createM1FromM0andFragAnnotation(trees[tree][specID],masses,sub_forms,carbon_pos,bounds,ppmError)[0]
                else:
                    spec = trees[tree][specID]
                returnDict[(id,mz,name,tree[4],tree[3])][specID] = collapseAsNeeded(spec,res)


        return returnDict

    """
    take spectra from m/z cloud after being read in by json and convert to dictionary
    """
    @staticmethod
    def convertSpectra2Vector(spectra):
        mz = spectra['IsolationWidth']["XPos"]
        spectrum = {x["MZ"]:0 for x in spectra["Peaks"]}
        for x in spectra["Peaks"]:
            spectrum[x["MZ"]] += x["Abundance"]
        return spectrum,mz

    @staticmethod
    def getTrees(trees, keys, calibration="recalibrated",library="reference"):
        output = {}
        for tree in trees:
            url = keys.SPECTRAURL.replace("TREENUMBER", str(tree[0]))
            url = url.replace("LIBRARY", library)
            querystring = {"stage": "2", "processing": calibration, "peaks": "true"}
            payload = ""
            response = 1
            while (response == 1):
                try:
                    response = requests.get(url, data=payload, headers=keys.HEADER, params=querystring, timeout=30)
                    if type(response) == type(
                            None) or "An error has occurred" in response.text or "Service Unavailable" in response.text or "unavailable" in response.text:
                        response = 1
                    else:
                        specTree = mzCloudPy.getAllSpectraInTree(
                            mzCloudPy.reformatSpectraDictList(json.loads(response.text)))
                        specTreeClean = {}
                        for spec in specTree:
                            if 1e6 * abs(specTree[spec][1]-tree[5])/tree[5] < 100:
                                specTreeClean[spec] = specTree[spec][0]
                        if len(specTreeClean) > 0:
                            output[tree] = specTreeClean
                except:
                   pass

        return output
    @staticmethod
    def reformatSpectraDictList(dictList):
        dat = {x["Id"]: {key: val for key, val in x.items() if key != "Id"} for x in dictList}
        return dat
    @staticmethod
    def getAllSpectraInTree(dat):
        return {id: mzCloudPy.convertSpectra2Vector(dat[id]) for id in dat}

    @staticmethod
    def getCompoundList(page,keys, pageSize=100, library="reference"):
        url = keys.COMPOUNDURL
        url = url.replace("LIBRARY", library)
        querystring = {"page": str(page), "pageSize": str(pageSize), "newer": "2000-09-03"}
        payload = ""

        response = 1
        while (response == 1):
            try:
                response = requests.get(url, data=payload, headers=keys.HEADER, params=querystring, timeout=timeout)
                if type(response) == type(
                        None) or "An error has occurred" in response.text or "Service Unavailable" in response.text:
                    response = 1
                else:
                    dat = mzCloudPy.reformatSpectraDictList(json.loads(response.text)["Items"])
            except:
                pass
        return dat

    @staticmethod
    def generateCompoundID2SpectralIDIndexedByM_ZStrict(numPerPage=100, key="none",library="reference"):
        # get number of items
        keys = Keys(key)
        totalCompounds = json.loads(requests.get(keys.COMPOUNDURL.replace("LIBRARY", library), data="",
                                                                 headers=keys.HEADER, params={"page": str(1),
                                                                                              "pageSize": "5",
                                                                                              "newer": "2000-09-03"}).text)["Total"]

        print(totalCompounds)
        numPages = int(np.ceil(float(totalCompounds) / numPerPage))
        linkage = {"Positive": [], "Negative": []}
        for page in range(1, numPages + 1):
            compounds = mzCloudPy.getCompoundList(page, keys,numPerPage, library=library)
            treeDict = {}
            for comp in compounds:
                try:
                    formula = compounds[comp]["InChI"].split("/")[1]
                    tmp = molmass.Formula(formula)
                    mz = tmp.isotope.mass
                    treeDict.update({(key, compounds[comp]["InChIKey"], compounds[comp]["SearchCompoundName"],formula,-1,mz): val for key, val in
                                 mzCloudPy.reformatSpectraDictList(compounds[comp]["SpectralTrees"]).items()})
                except:
                    pass
            for tree in treeDict:
                try:
                    polarity = treeDict[tree]["Polarity"]
                    if polarity == "Positive":
                        mz = tree[-1] + 1.0073
                    elif polarity == "Negative":
                        mz = tree[-1] - 1.0073
                    else:
                        print("Error",polarity)
                    tree = tuple(list(tree[:-1]) + [mz])
                    linkage[polarity].append(tree)
                except:
                    pass

            print(float(page) * numPerPage / totalCompounds)
        pkl.dump(linkage, open("../src/DecoID/mzCloudCompound2TreeLinkage_InChI" + library + ".pkl", "wb"), pkl.HIGHEST_PROTOCOL)







