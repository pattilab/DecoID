import customDBpy
import mzCloudPy
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
sys.path.append(os.path.join(application_path,"..","src/"))
from multiprocessing import Queue, Process,Manager,Value,Lock,set_start_method
from threading import Thread
import time
import gzip
import pickle as pkl
import dill
import pandas as pd
from MS2Search import *
import uuid

class DecoID():
    def __init__(self,useCache,libFile,useAuto = False,numCores=1,resolution = 2,label=""):
        self.useCache = useCache
        self.libFile = libFile
        self.useDBLock = False
        if ".tsv" in libFile:
            try:
                f = open(libFile,"r",encoding = "utf-8")
                data = [x.replace(",","_").rstrip().split("\t") for x in f.readlines()]
                header = data[0]
                data = [{k:replaceAllCommas(val) for k,val in zip(header,x)} for x in data[1:]]
                for x in data:
                    x["m/z"] = float(x["m/z"])
                    x["spectrum"] = createDictFromString(x["spectrum"],resolution)

                posData = {x["id"]:x for x in data if "ositive" in x["mode"]}
                negData = {x["id"]:x for x in data if "egative" in x["mode"]}
                self.library = {"Positive":posData,"Negative":negData}
            except:
                print("bad library file: ", libFile)
                return -1
            self.lib = customDBpy
            self.useCache = False
            self.cachedReq = "none"
            self.ms_ms_library = "custom"
            self.useAuto = useAuto
            pkl.dump(self.library, open(libFile.replace(".tsv", ".db"),"wb"))
            print("Library loaded successfully: " + str(
            len(self.library["Positive"]) + len(self.library["Negative"])) + " spectra found")

        elif ".msp" in libFile:
            try:
                mz = -1
                id = -1
                cpdid = -1
                polarity = -1
                exactMass = -1
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
                                if polarity != -1:
                                    if mz == -1:
                                        if polarity == "Negative":
                                            mz = exactMass - 1.0073
                                        else:
                                            mz = exactMass + 1.0073
                                    if cpdid == -1:
                                        cpdid = name
                                    if id != -1 and len(spectrum) > 0:
                                        self.library[polarity][id] = {"id":id,"cpdID":cpdid,"name":replaceAllCommas(name),"mode":polarity,"spectrum":spectrum,"m/z":np.round(mz,4)}
                                        good = True


                        else:
                            first = False
                        name = line.split(": ")[1]
                        mz = -1
                        id = -1
                        cpdid = -1
                        polarity = -1
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
                    if "PrecursorMZ: " in line:
                        mz = float(line.split(": ")[1])
                    if "Ion_mode: " in line:
                        temp = line.split(": ")[1]
                        if temp == "P":
                            polarity = "Positive"
                        elif temp == "N":
                            polarity = "Negative"
                    if looking and len(line.split()) == 2:
                        spectrum[float(line.split()[0])] = float(line.split()[1])
                    if "Num Peaks" in line:
                        looking = True
                print("Library loaded successfully: " + str(len(self.library["Positive"]) + len(self.library["Negative"])) + " spectra found")
                pkl.dump(self.library,open(libFile.replace(".msp",".db"),"wb"))
            except:
               print(sys.exc_info())
               print("bad library file: ", libFile)
               return -1
            self.lib = customDBpy
            self.useCache = False
            self.cachedReq = "none"
            self.ms_ms_library = "custom"
            self.useAuto = useAuto

        elif ".db" in libFile:
            self.lib = customDBpy
            self.useCache = False
            self.cachedReq = "none"
            self.ms_ms_library = "custom"
            self.useAuto = useAuto
            self.library = pkl.load(open(libFile,"rb"))
            print("Library loaded successfully: " + str(
                len(self.library["Positive"]) + len(self.library["Negative"])) + " spectra found")

        else:
            self.useDBLock = True
            self.lib = mzCloudPy
            self.ms_ms_library = "mzCloud"
            if useCache:
                manager = Manager()
                self.cachedReq = manager.dict()
                good = True
                try:
                    fp = gzip.open(self.lib.CACHEFILE + ".cache", "rb")
                    cacheDict = pkl.load(fp)
                except:
                    good = False
                if good:
                    chunkSize = 500
                    allKeys = list(cacheDict.keys())
                    if len(allKeys) > 0:
                        for smallChunks in splitList(allKeys, max([1, int(len(allKeys) / chunkSize)])):
                            self.cachedReq.update({key: cacheDict[key] for key in smallChunks})
                    print(len(self.cachedReq), " cached trees")
                    fp.close()
            else:
                self.cachedReq = "none"
            self.library = ["reference"]
            self.useAuto = useAuto


            if useAuto:
                self.library.append("autoprocessing")
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


    def readData(self,filename,resolution,peaks,DDA,massAcc,offset=.65,peakDefinitions = "",tic_cutoff=0):
        samples, ms1 = readRawDataFile(filename, MAXMASS, resolution, peaks,offset=offset,tic_cutoff=tic_cutoff,ppmWidth=massAcc)
        if samples != -1:
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
            print(len(self.samples), " MS2 spectra detected")
            fileending = "." + filename.split(".")[-1]
            self.filename = filename.replace(fileending, "")
            self.filename = self.filename.replace('"', "")
            self.clusters = {}
            self.data2WriteFinal = []
            self.outputDataFinal = {}

            if peakDefinitions != "":
                self.peakData = pd.read_csv(peakDefinitions)
                self.peakData = self.peakData[["mz","rt_start","rt_end"]]
                dfIndex = list(self.peakData.index.values)

                if self.DDA:
                    i = -1
                    goodSamps = []
                    for samp in samples:
                        i += 1
                        for index in dfIndex:#,row in self.peakData.iterrows():
                            if abs(samp["center m/z"]-self.peakData.at[index,"mz"])/self.peakData.at[index,"mz"]*1e6 < self.massAcc:
                                if samp["rt"] >= self.peakData.at[index,"rt_start"] and samp["rt"] <= self.peakData.at[index,"rt_end"]:
                                    samp["group"] = index
                                    goodSamps.append(i)
                                    break
                    self.samples = [self.samples[x] for x in goodSamps]
                    numGroups = len(set([samp["group"] for samp in self.samples]))
                    print("Number of compounds with acquired MS2: ",numGroups)
                    print("Number of spectra to deconvolve: ",len(self.samples))

                else:
                    samplesDict = {x["id"]: x for x in samples}
                    mzGroups = clusterSpectraByMzs([x["center m/z"] for x in samples], [x["id"] for x in samples],
                                                   [x["rt"] for x in samples],list(range(len(samples))) ,self.massAcc)
                    index = 0
                    for mz in mzGroups:
                        for m, id in mzGroups[mz]:
                            samplesDict[id]["group"] = index
                        index += 1

            else:
                index = 0
                for samp in samples:
                    samp["group"] = index
                    index += 1
                self.peakData = pd.DataFrame.from_dict({i:{"mz":samp["center m/z"],"rt_start":samp["rt"]-.05,"rt_end":samp["rt"]+.05} for i,samp in zip(range(len(samples)),samples)},orient="index")

    def readMS_DIAL_data(self,file,mode,massAcc,peakDataFile):
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

    def identifyUnknowns(self,resPenalty=100,percentPeaks=0.01,iso=False,ppmThresh = 10,dpThresh = 20):
        try: self.samples
        except: print("datafile not loaded");return -1
        self.iso = iso
        self.percentPeaks = percentPeaks

        self.resPenalty = resPenalty

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
                if (x["percentContamination"] < .2) and x["signal"] > sigThresh:
                    samplesToGo.append(x)
                else:
                    samplesLeft.append(x)

            data2Write = []
            outputData = {}
            t = Thread(target=self.recursiveRunQ, args=(q,data2Write,outputData))
            t.start()

            self.runSamples(samplesToGo,q)

            t.join()


            unknownGroupIDs = []
            unknownThreshold = dpThresh
            unknownPPMThreshold = ppmThresh
            for result, centerMz, id, rt, s2n, indices, numComponents, fragments, decoSpec,components in data2Write:
                if not any(result[x][0] >= unknownThreshold and ((result[x][4] - x[3]) * (10 ** 6)) / x[
                    3] < unknownPPMThreshold for x in result):
                    unknownGroupIDs.append(id)

            unknownSamples = []
            for id in unknownGroupIDs:
                for samp in samplesToGo:
                    if samp["group"] == id:
                        unknownSamples.append(samp)


            self.clusters = clusterSpectra(unknownSamples,self.peakData)

            print("Number of Predicted Unknown Compounds: ", len(self.clusters))



    def searchSpectra(self,verbose,resPenalty = 100,percentPeaks=0.01,iso=False,threshold = 0.0):
        try: self.samples
        except: print("datafile not loaded");return -1
        q = Queue(maxsize=1000)

        self.iso = iso
        self.percentPeaks = percentPeaks

        self.resPenalty = resPenalty
        if len(self.ms1) < 1:
            self.recursive = False
            self.peaks = False
            # if resPenalty > 10:
            #     self.resPenalty = 10

        if not self.peaks and not self.DDA:
            self.DDA = True


        if type(self.samples) != type(-1):
            self.paramSuffix = ""
            for x in [self.useAuto, self.recursive, self.iso, self.peaks]:
                self.paramSuffix += "_" + str(int(x))
            success = False
            error = True
            while (not success):
                #try:
                    outfile = open(self.filename+ self.label + "_decoID" + ".csv", "w")
                    success = True
                # except:
                #     if error:
                #         print("Error: please close " + self.filename + self.label + "_decoID.csv")
                #         error = False
                #         time.sleep(2)

            outfile.write(
                "#scanID,isolation_center_m/z,rt,compound_m/z,DB_Compound_ID,Compound_Name,DB_Spectrum_ID,dot_product,ppm_Error,Abundance,ComponentID\n")
            t = Thread(target=self.runQ, args=(
            q,outfile,threshold,self.outputDataFinal,self.filename + self.label + ".DecoID",
            self.filename + self.label + "_scanInfo.csv", verbose))
            t.start()

            t2 = Thread(target=self.sendCompleteSamples, args=(self.data2WriteFinal,q))
            t2.start()

            self.runSamples(self.samples,q)


            t.join()
            t2.join()
        if self.useCache:
            chunkSize = 500
            fp = gzip.open(self.lib.CACHEFILE + ".cache", "wb")
            allKeys = list(self.cachedReq.keys())
            newDict = dict()
            for smallChunks in splitList(allKeys, max([1, int(len(allKeys) / chunkSize)])):
                newDict.update({key: self.cachedReq[key] for key in smallChunks})
            pkl.dump(newDict, fp)
            fp.close()

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
        outputScanFile = open(scanInfoFileName, "w")
        outputScanFile.write("#scanID,Signal to Noise Ratio,numComponents,componentID,componentAbundance,componentMz,rt,spectrum\n")
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
                        outputScanFile.write(str(id) + "," + str(s2n) + "," + str(numComponents) + "," + str(c[0]) + "," +str(c[2]) + "," + str(c[1]) + "," + str(rt) + ",")
                        for ind,i in zip(indices,components[c]):
                            if i > 0:
                                outputScanFile.write(str(indexRef[ind]) + ":" + str(i) + " ")
                        outputScanFile.write("\n")
                    update = False
                    if id not in outputDataFile:
                        outputDataFile[id] = {}
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

                            try:
                                [outfile.write(delimiter + z) for z in
                                 [str(rt), str(x[3]), str(x[1]), x[0], str(x[2]), str(result[x][0]),
                                  str((x[3] - tempMz) * 1e6 / tempMz), str(x[4]),str(componentName)]]
                            except:
                                [outfile.write(delimiter + z) for z in
                                 [str(x[1]), str(x[2]), str(result[x][0]),
                                  str((x[3] - tempMz) * 10 ** 6 / tempMz), str(x[4]),str(componentName)]]
                            outfile.write("\n")
                            if filename != "None" and update:
                                outputDataFile[id]["Hits"][x] = result[x]
                index += 1
                if type(verbose) == type(str()):
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
        # if self.recursive:
        #     tempObject = DecoID.from_DecoID(self)
        #     tempObject.filename+="_unknowns_"
        #     tempObject.samples = []
        #     DecoID.writeObj(tempObject,tempObject.filename+".dill")
        groups = list(set([x["group"] for x in self.samples]))
        numGroupsPerFile = int(len(groups)/numberOfFiles)
        self.data2WriteFinal = []
        self.outputDataFinal = {}
        groupBegin = 0
        processes = []
        groupEnd = groupBegin + numGroupsPerFile
        for num in range(numberOfFiles):
            tempObject = DecoID.from_DecoID(self)
            tempObject.filename += "_" + str(num) + "_"
            if num != numberOfFiles - 1:
                group = groups[groupBegin:groupEnd]
            else:
                group = groups[groupBegin:]
            groupBegin = groupEnd
            groupEnd = groupBegin + numGroupsPerFile
            tempObject.samples = [samp for samp in tempObject.samples if samp["group"] in group]
            DecoID.writeObj(tempObject,tempObject.filename+".dill")

    @staticmethod
    def combineResults(filenames,newFilename,endings = ["_scanInfo.csv","_decoID.csv",".DecoID"]):
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
                # else:
                #     print(end)
                #     fh = gzip.open(goodFiles[0]+end, "rb")
                #     [samples,ms1,peakData,baseDict] = pkl.load(fh)
                #     fh.close()
                #     if len(goodFiles) > 1:
                #         for f in goodFiles[1:]:
                #             fh = gzip.open(f + end, "rb")
                #             [samp,_,_,tempDict] = pkl.load(fh)
                #             baseDict.update(tempDict)
                #             samples += samp
                #             fh.close()
                #     fh = gzip.open(newFilename+end,"wb")
                #     pkl.dump([samples,ms1,peakData,baseDict],fh)
        else:
            print("no files detected")
    def runSamples(self,samples,q,numConcurrentMzGroups=20):

        numProcesses = Value('i', 0)
        lock = Lock()
        processes = []
        dbLock = Lock()
        availableToGrab = Value('i',1)
        samplesDict = {x["id"]:i for x,i in zip(samples,range(len(samples)))}
        keys = list(samplesDict.keys())
        mzGroups  = clusterSpectraByMzs([samples[samplesDict[k]]["center m/z"] for k in keys],[samples[samplesDict[k]]["id"] for k in keys],[samples[samplesDict[k]]["rt"] for k in keys],[samples[samplesDict[k]]["group"] for k in keys],self.massAcc)
        mzGroupsWSamples = {x:{"samples":[]} for x in mzGroups}
        for mz in mzGroups:
            for m,id in mzGroups[mz]:
                mzGroupsWSamples[mz]["samples"].append(samples[samplesDict[id]])

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

        def startProc(numP, l,samples,trees,possCompounds,possIsotopes):
            spectra = [samp["spectra"] for samp in samples]
            toAdd = {key: value for key, value in clusters.items() if
                     any(key[1] < samp["rt"] and key[2] > samp["rt"] for samp in samples)}
            p = Process(target=DecoID.processGroup,
                        args=(
                            samples,group,trees, spectra, possCompounds, possIsotopes, self.iso,self.DDA,self.massAcc,self.peaks,q,self.resPenalty,self.resolution,toAdd,self.peakData,lowerBound,upperBound))
            #t = Thread(target=startProc, args=(p, numProcesses, lock))
            while not checkRoom(numProcesses, lock):
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
                trees, possCompounds, possIsotopes = self.lib.getCanidateCompoundTrees(mode, upperBound, lowerBound, self.iso,
                                                                          self.library,self.cachedReq)
        else:
            trees, possCompounds, possIsotopes = self.lib.getCanidateCompoundTrees(mode, upperBound, lowerBound,
                                                                                   self.iso,
                                                                                   self.library, self.cachedReq)
        featureGroups = {x["group"]:[] for x in allSamples}
        [featureGroups[x["group"]].append(x) for x in allSamples]



        threads = []
        for group,samples in featureGroups.items():
            t = Thread(target=startProc,args=(numProcesses,lock,samples,trees,possCompounds,possIsotopes))
            t.start()
            threads.append(t)

        for t in threads:
            t.join()

    @staticmethod
    def processGroup(samples,group,trees, spectra, possCompounds, possIsotopes, iso,DDA,massAcc,peaks,q,resPenalty,resolution,toAdd,peakData,lowerbound,upperbound):

        [metIDs, spectraTrees, spectraIDs, matrix, masses, metNames, indicesAll,
         reduceSpec] = getMatricesForGroup(trees, spectra, possCompounds, possIsotopes, iso,
                                           # make get matrices for group independent of library source
                                           resolution, toAdd)
        isoIndices = [x for x in range(len(metIDs)) if type(metIDs[x]) == type(tuple())]

        results = []
        samples2Go = list(samples)
        reduceSpec2Go = list(reduceSpec)

        for sample, spectrum in zip(samples2Go, reduceSpec2Go):
            results.append(DecoID.processSample(sample, spectrum, masses, matrix, massAcc, peaks,resPenalty,isoIndices))

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

            scores, components = scoreDeconvolution(combinedSpectrum, matrix, resVector, metIDs, spectraTrees,
                                                    masses, centerMz, DDA, massAcc)

            result = {(met, specID, mass, name, safeDivide(comp, sum(resVector))): coeff for
                      met, name, specID, coeff, mass, comp in
                      zip(metIDs, metNames, spectraIDs, scores, masses, resVector)}

            # combine isotopes together
            result = cleanIsotopes(result)

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
            rts = [samp["rt"] for samp in samples]
            for index,row in peakData.iterrows():
                if row["mz"] >= lowerbound and row["mz"] <= upperbound:
                    goodIndices = [x for x in range(len(rts)) if rts[x] >= row["rt_start"] and rts[x] <= row["rt_end"]]
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

                        scores, components = scoreDeconvolution(combinedSpectrum, matrix, resVector, metIDs, spectraTrees,
                                                                masses, centerMz, DDA, massAcc)

                        result = {(met, specID, mass, name, safeDivide(comp, sum(resVector))): coeff for
                                  met, name, specID, coeff, mass, comp in
                                  zip(metIDs, metNames, spectraIDs, scores, masses, resVector)}

                        # combine isotopes together
                        result = cleanIsotopes(result)

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

            #get peaks in window
            #parse out by retention time where peaks occur
            #repeat the above with the targeted m/z
            #no peaks => no results

    @staticmethod
    def processSample(sample,spectrum,masses,matrix,massAcc,peaks,resPenalty,isoIndices):

        centerMz = sample["center m/z"]
        lowerBound = sample["lower m/z"]  # lowest m/z in isolation window
        upperBound = sample["higher m/z"]  # highest

        def checkScan(mass,index):
            if index in isoIndices:
                return inScanIso(sample["fragments"],mass,massAcc)
            else:
                return inScan(sample["fragments"],mass,massAcc)

        if len(matrix) > 0:
            goodIndices = [x for x in range(len(masses)) if masses[x] < upperBound and masses[x] > lowerBound and ((not peaks) or checkScan(masses[x],x))]
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
            #q.put([res,s2n,decoSpec])
            return [res,s2n,decoSpec]
        else:
            #q.put([[],0,[0 for _ in spectrum]])
            return [[],0,[0 for _ in spectrum]]


            # score the components found in the deconvolution

        #    numComponents = len([x for x in res if x > 0])

        #
        #     scores,components = metaboliteSensitivity(spectra[-1], matrices[-1], res, metIDs[-1], spectraTrees[-1],masses[-1],centerMz,DDA,massAcc,fragments)
        #
        #     result = {(met, specID, mass, name, safeDivide(comp, sum(res))): coeff for
        #               met, name, specID, coeff, mass, comp in
        #               zip(metIDs[-1], metNames[-1], spectraIDs[-1], scores, masses[-1], res)}
        #
        #     # combine isotopes together
        #     result = cleanIsotopes(result)
        #
        #     # output result
        #     success = False
        #     while not success:
        #         try:
        #             q.put([result, centerMz, id, rt, s2n, indices[-1], numComponents, fragments, decoSpec,components ], timeout=5)
        #             success = True
        #         except:
        #             pass
        # else:
        #     q.put([[], centerMz, id, rt, 0, [], 0, [], [0 for _ in spectra[0]],{}])


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









