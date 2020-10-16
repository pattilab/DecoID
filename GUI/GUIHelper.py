from tkinter import *
from tkinter import font
from tkinter import filedialog,ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
from matplotlib.figure import Figure
import time
import gzip
from DecoID.DecoID import *
from multiprocessing import Queue,Process
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.dpi'] = 300
import os
import numpy as np
import pickle as pkl

LARGE_FONT = ("Verdana", 12)
RESOLUTIONS = [-1,0,1,2]
METHODS = ["dot","eV","NCE"]
PPMFILTER = [False,True]
RT = [False,True]
MZCLOUDSELECTION = ["reference","autoprocessing","none","none","none"]

if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
    DATABASESELCTION = ["none","none", os.path.join(application_path,"MoNA-export-LC-MS-MS_Spectra.db"), os.path.join(application_path,"HMDB_experimental.db"),
                        "custom"]
elif __file__:
    application_path = os.path.dirname(__file__)
    DATABASESELCTION = ["none","none", os.path.join(application_path,"../databases/MoNA-export-LC-MS-MS_Spectra.db"), os.path.join(application_path,"../databases/HMDB_experimental.db"),
                        "custom"]

class decoIDSearch(Tk):

    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        Tk.wm_title(self, "decoID Search")

        container = Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}

        startFrame = StartPage(container,self)
        #startFrame.grid(row=0, column=0, sticky="nsew")

        # game = Game(container,self,initializeSara,otherPlayers)
        # game.grid(row=0, column=0, sticky="nsew")

        # self.frames[Game] = game
        self.frames[StartPage] = startFrame

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()



class Checkbar(Frame):
   def __init__(self, parent=None, picks=[], side=LEFT, anchor=W,label = ""):
      Frame.__init__(self, parent)
      self.vars = []
      label = Label(self, text=label)
      label.pack()
      for pick in picks:
         var = IntVar()
         chk = Checkbutton(self, text=pick, variable=var)
         chk.pack(side=side, anchor=anchor, expand=YES)
         self.vars.append(var)
   def state(self):
      return map((lambda var: var.get()), self.vars)

class Radiobar(Frame):
   def __init__(self, parent=None, picks=[], side=LEFT, anchor=W,label = ""):
      Frame.__init__(self, parent)
      label = Label(self, text=label)
      label.grid(row=0,column=0,padx=10)
      self.var = IntVar()
      for pick,num in zip(picks,range(len(picks))):
         var = IntVar()
         chk = Radiobutton(self, text=pick, variable=self.var,value=num)
         chk.grid(row=0, column=1+num,padx=5)
   def state(self):
      return self.var.get()

def browseForFile(filename,types):
    filename.set(filedialog.askopenfilename(filetypes=types))

def performSearch(filename,numCores,recursive,iso,peaks,dtype,massAcc,custOrMzCloud,libFile,searchMethod,filenamePeak,offset,useRT,rtTol,key,fragInt):
    file = filename.get()
    libraryFile = libFile.get()
    peakFile = filenamePeak.get()
    if peakFile == "No File Selected":
        peakFile = ""

    fragIntThresh = float(fragInt.get())

    try:
        apikey = open(key.get()).readline().rstrip()
    except:
        apikey = "none"

    db = DATABASESELCTION[custOrMzCloud.state()]
    if db != "custom":
        libraryFile = db

    mzCloudLib = MZCLOUDSELECTION[custOrMzCloud.state()]


    usert = RT[useRT.state()]
    if useRT:
        rtTolerance = float(rtTol.get())

    else:
        rtTolerance = float("inf")

    threshold = 0#scale.get()

    useIso = 0
    usePeaks = 0
    doDeco = 1
    DDA = 1

    if not PPMFILTER[searchMethod.state()]:
        doDeco = 0
    #if not PPMFILTER[library.state()]:
    #    useAuto = 1
    useRec = 0
    if not PPMFILTER[recursive.state()]:
        useRec = 1
    if not PPMFILTER[iso.state()]:
        useIso = 1
    if not PPMFILTER[peaks.state()]:
        usePeaks = 1
    if PPMFILTER[dtype.state()]:
        DDA = 0
    MA = massAcc.get()

    fragThresh = 0#scaleFragment.get()


    thr = threading.Thread(target=createWaitingBox,args=(file,threshold,fragThresh,numCores.get(),mzCloudLib,useRec,useIso,usePeaks,DDA,MA,libraryFile,doDeco,peakFile,offset,rtTolerance,apikey,fragIntThresh))
    thr.start()


def createWaitingBox(file,threshold,fragThresh,numCores,mzcloudLib,useRec,useIso,usePeaks,DDA,massAcc,libFile,doDeco,peakFile,offset,rtTolerance,key,fragIntThresh):
    def cancelSearch(proc,var,top):
        proc.kill()
        var.set("Search Canceled")
        top.destroy()
    top = Toplevel()
    var = StringVar()
    var.set("Searching... \n\n\n")
    lab = Label(top,textvariable=var, font=LARGE_FONT)
    lab.pack()
    progress = ttk.Progressbar(top, orient=HORIZONTAL, length=500, mode='determinate')
    progress.pack()
    Label(top,text="\n\n",font=LARGE_FONT).pack()

    q = Queue()
    progress['value'] = 1.0
    numCores = int(numCores)
    def runDeco():
        decID = DecoID(libFile, mzcloudLib, numCores,api_key=key)
        decID.readData(file ,2,usePeaks,DDA,massAcc,peakDefinitions=peakFile,offset=offset,frag_cutoff=fragIntThresh)

        if doDeco:
            if DDA: lam = 1
            else: lam = 100
            if useRec and DDA:
                decID.identifyUnknowns(iso=useIso)
            decID.searchSpectra(q,lam,fragThresh,useIso,threshold,rtTol=rtTolerance)
        else:
            decID.searchSpectra(q,np.inf,fragThresh,useIso,threshold,rtTol=rtTolerance)

    p = Process(target=runDeco(),args=())


    p.start()
    # popen = subprocess.Popen(
    #     cmd, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1)

    cancelButton = ttk.Button(top, text="Cancel",
                              command=lambda : cancelSearch(p,var,top))
    cancelButton.pack()


    while( p.is_alive() or not q.empty()):#and for path in execute(popen):
        if not q.empty():
            q.get()
            progress["value"] = min([100,2+progress["value"]])
        else:
            time.sleep(5)

    var.set("Search Complete: output in raw datafile directory")

def displaySpectraWindow(filename):
    resolution = 1
    result = readRawDataFile('"' + filename.get() + '"', maxMass=3000, resolution=[resolution])
    displaySpectra(result.__iter__())


def plotSpectra(minMass,maxMass,res,spectras,cs=["blue","red"],labels=["Original","Deconvolved"],title="",specificPlot = -1):
    plt.figure()
    if type(spectras[0]) != type(list()):
        spectras = [spectras]
    spectras = [normalizeSpectra(spec,"max") for spec in spectras]
    cs = cs[:len(spectras)]
    coefficient = 1
    mzs = np.linspace(minMass,maxMass,int(maxMass/res))
    # smallestmz = mzs[min([min([x for x in range(len(spectra)) if spectra[x] >= 1e-3]) for spectra in spectras])]
    # largestmz = mzs[max([max([x for x in range(len(spectra)) if spectra[x] >= 1e-3]) for spectra in spectras])]
    goodInd = flatten([[x for x in range(len(spectra)) if spectra[x] >= 1e-3] for spectra in spectras])
    if len(goodInd) > 0:
        smallestmz = mzs[min(goodInd)]
        largestmz = mzs[max(goodInd)]
    else:
        smallestmz = minMass
        largestmz = maxMass

    for spectra,c,label in zip(spectras,cs,labels):
        for m,i in zip(mzs,spectra):
            if i >= 1e-3:
                if specificPlot == -1:
                    plt.plot([m,m],[0,i*coefficient],c=c)
                else:
                    specificPlot.plot([m,m],[0,i*coefficient],c=c)
        if specificPlot == -1:
            plt.plot([largestmz+1,largestmz+1],[0,0],c=c,label=label)
        else:
            specificPlot.plot([largestmz+1,largestmz+1],[0,0],c=c,label=label)
        coefficient *= -1
    if specificPlot == -1:
        plt.plot([smallestmz,largestmz+1],[0,0],c="black",linewidth=.5)
    else:
        specificPlot.plot([smallestmz, largestmz+1], [0, 0], c="black", linewidth=.5)
    maxIntenstiy = max(flatten(spectras))
    interval = maxIntenstiy/2
    ticks = []
    current = 0.0
    ticks.append(current)
    while(current < maxIntenstiy):
        current += interval
        ticks.append(np.round(current,4))
    ticksReverse = ticks[1:]
    ticksReverse.reverse()


    if specificPlot == -1:
        plt.yticks([-1*x for x in ticksReverse]+ticks,ticksReverse+ticks)
        plt.title(title)
        plt.xlabel("m/z")
        plt.ylabel("Relative Intensity")
        plt.legend()
        if len(spectras) == 1:
            plt.ylim((0,1))
    else:
        specificPlot.set_yticks([-1*x for x in ticksReverse]+ticks)
        specificPlot.set_yticklabels([str(x) for x in ticksReverse+ticks])
        specificPlot.set_title(title)
        specificPlot.set_xlabel("m/z")
        specificPlot.set_ylabel("Relative Intensity")
        specificPlot.legend()
        if len(spectras) == 1:
            specificPlot.set_ylim((0,1))




def displaySpectra(result):
    sample = next(result)
    top = Toplevel()
    w = 1000
    h = 1400
    x = 400
    y = 100
    top.geometry('%dx%d+%d+%d' % (w, h, x, y))
    f = Figure(figsize=(10, 10), dpi=100)
    a = f.add_subplot(111)
    spectra = np.zeros(3000).tolist()
    for mz,i in sample["spectra"][0].items():
        spectra[int(np.round(mz,0))] += i

    plotSpectra(0, maxMass=MAXMASS, res=10 ** (-1 * 1), spectras=[spectra], specificPlot=a,
                title="m/z = " + str(sample["center m/z"]) + str("   rt = " + str(sample["rt"])))
    canvas = FigureCanvasTkAgg(f, top)
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
    nextFigureButton = ttk.Button(top, text="Next",
                                 command=lambda: closeAndRestart(top,lambda:displaySpectra(result)))
    nextFigureButton.pack()

def closeAndRestart(pane,command):
    pane.destroy()
    try:
        command()
    except:
        pass


def visualizeResultsStart(self,filename,scoreThresh,ppmThresh,usePPM,scanNumber,recursive,iso,peaks,DDA):

    top = Toplevel()
    if PPMFILTER[usePPM]:
        ppmThresh = np.inf

    if PPMFILTER[DDA]:
        DDA = False
    else:
        DDA = True
    Label(top, text="Clustering...", font=LARGE_FONT).pack()

    thr = threading.Thread(target=visualizeResults,
                           args=(self,filename,top,scoreThresh,ppmThresh,scanNumber,DDA))
    thr.start()

def visualizeResults(self, filename, oldFrame, scoreThresh, ppmThresh, groupNumber, DDA):
    name = filename.get()
    ending = name.split(".")[-1]
    name = name.replace("." + ending,".DecoID")

    try:
        [samples,peakData,ms1,results] = pkl.load(gzip.open(name,"rb"))
        #peaks = bool(int(prefix.split("_")[-1]))
        #samples,ms1 = readRawDataFile('"' + filename.get() + '"', maxMass=3000, resolution=[2], useMS1=peaks)
        if len(ms1) < 1: peaks = False
        if DDA:
            isoWidth = abs(samples[0]["lower m/z"] - samples[0]["higher m/z"])
            if groupNumber != "":
                groupNumber = int(groupNumber)
                samples = [x for x in samples if str(x["group"]) == str(groupNumber)]
            # spectra = {x["id"]: x["spectra"] for x in samples}
            # if peaks:
            #     ms1 = {x["id"]: x["ms1"] for x in samples}
            # else:
            #     ms1 = {x["id"]:[] for x in samples}
            clusteredSamples = clusterSpectra(samples,peakData)
            rts = [c["rtWindow"] for c in clusteredSamples]
            mzC = [c["m/z"] for c in clusteredSamples]
            ordering = list(range(len(mzC)))
            ordering.sort(key=lambda x:mzC[x])
            mzC = [mzC[x] for x in ordering]
            rts = [rts[x] for x in ordering]
            clusteredSamples = [clusteredSamples[x] for x in ordering]
            isoBounds = [[mz-isoWidth/2,mz+isoWidth/2] for mz in mzC]
        else:
            isoWidth = np.mean([np.abs(s["higher m/z"] - s["lower m/z"]) for s in samples])
            if groupNumber != "":
                groupNumber = int(groupNumber)
                s = [x for x in samples if x["lower m/z"] <= peakData.at[groupNumber,"mz"] and x["higher m/z"] >= peakData.at[groupNumber,"mz"] and x["rt"] >= peakData.at[groupNumber,"rt_start"] and x["rt"] <= peakData.at[groupNumber,"rt_end"]]
                if len(s) > 0:
                    rts = [[peakData.at[groupNumber,"rt_start"],peakData.at[groupNumber,"rt_end"]]]
                    mzC = [peakData.at[groupNumber,"mz"]]
                    spec = s[0]["spectra"]
                    if len(s) > 1:
                        for spec2 in s[1:]:
                            spec = mergeSpectrum(spec, spec2["spectra"])
                    clusteredSamples = [
                        {"group": groupNumber, "spectrum": spec, "m/z": np.mean([x["center m/z"] for x in s]),
                         "rtWindow": rts[-1]}]
                isoBounds = [[s[0]["lower m/z"],s[0]["higher m/z"]]]
            else:
                rts = []
                mzC = []
                isoBounds = []
                clusteredSamples = []
                for groupNumber,row in peakData.iterrows():
                    s = [x for x in samples if
                               x["lower m/z"] <= peakData.at[groupNumber, "mz"] and x["higher m/z"] >= peakData.at[
                                   groupNumber, "mz"] and x["rt"] >= peakData.at[groupNumber, "rt_start"] and x[
                                   "rt"] <= peakData.at[groupNumber, "rt_end"]]
                    if len(s) > 0:
                        spec = s[0]["spectra"]
                        if len(s) > 1:
                            for spec2 in s[1:]:
                                spec = mergeSpectrum(spec, spec2["spectra"])
                        rts.append([row["rt_start"],row["rt_end"]])
                        mzC.append(row["mz"])
                        clusteredSamples.append({"group": groupNumber, "spectrum": spec, "m/z": np.mean([x["center m/z"] for x in s]),
                         "rtWindow": rts[-1]})
                        isoBounds.append([s[0]["lower m/z"], s[0]["higher m/z"]])
                ordering = list(range(len(mzC)))
                ordering.sort(key=lambda x: mzC[x])
                mzC = [mzC[x] for x in ordering]
                rts = [rts[x] for x in ordering]
                isoBounds = [isoBounds[x] for x in ordering]
                clusteredSamples = [clusteredSamples[x] for x in ordering]
        featureLabels = [str(c["group"])+ ": m/z =  " + str(np.round(m,4)) + "| rt = ["+str(np.round(r[0],2))+","+str(np.round(r[1],2))+"]" for id,m,r,c in zip(range(len(mzC)),mzC,rts,clusteredSamples)]
        featureLabelDict = {key:val for key,val in zip(featureLabels,range(len(featureLabels)))}
        oldFrame.destroy()
        top = Toplevel()
        Label(top, text="Features Clustered", font=LARGE_FONT).pack()
        selection = StringVar()
        selection.set(featureLabels[0])
        featureLabels = featureLabels
        menu = OptionMenu(top,selection,*featureLabels)
        Label(top,text="Choose Feature").pack()
        menu.pack()

        displayButton = ttk.Button(top, text="Display Results",command = lambda :displayHitsForCluster(clusteredSamples[featureLabelDict[selection.get()]],results,selection.get(),filename.get(),ppmThresh,scoreThresh,DDA,ms1,isoBounds[featureLabelDict[selection.get()]]))
        displayButton.pack()


    except:
        top = Toplevel()
        Label(top, text="Datafile does not exist for selected file, please run the search first", font=LARGE_FONT).pack()


def displayHitsForCluster(cluster,results,feature,filename,ppmThresh,scoreThresh,DDA,ms1,isoBounds):

    relevant = {key:val for key,val in results.items() if key == cluster["group"]}
    newDict = {}
    ms1Deco = {}
    for id in relevant:
        comps = relevant[id]["components"]
        ms1Deco[id] = [[c[1],c[3]] for c in comps]
        ms1OI = [val for key,val in ms1.items() if key > cluster["rtWindow"][0] and key < cluster["rtWindow"][1]]
        if len(ms1OI) > 0:
            ms1Spec = {key:val for key,val in ms1OI[0].items() if key > isoBounds[0] and key < isoBounds[1]}
            for s in ms1OI[1:]:
                ms1Spec = mergeSpectrum(ms1Spec,{key:val for key,val in s.items() if key > isoBounds[0] and key < isoBounds[1]})
        else:
            ms1Spec = {}
        # if not DDA:
        #     tempMz = [f - float(relevant[id]["center m/z"]) for f in relevant[id]["fragments"]]
        #     tempMz.sort(key=lambda x: abs(x))
        #     tempMz = tempMz[0]
        # else:
        #     tempMz = relevant[id]["center m/z"]
        for hit in relevant[id]["Hits"]:
            tempMz = relevant[id]["center m/z"]
            if (abs(float(hit[2]) - float(tempMz))/float(tempMz)) * 10 **6 < ppmThresh and relevant[id]["Hits"][hit][0] > scoreThresh:
                newDict[(id,relevant[id]["rt"],tempMz,tuple(relevant[id]["index"]))+hit] = [relevant[id]["decoSpec"],ms1Spec] + relevant[id]["Hits"][hit]
    order = list(newDict.keys())
    order.sort(key=lambda x:newDict[x][2],reverse=True)
    result = [list(x) + newDict[x] for x in order]
    if len(result) > 0:
        displayHit(result.__iter__(),feature + " ("+str(len(result)) + " Hits)",cluster,filename,ms1Deco,isoBounds)
    else:
        top = Toplevel()
        Label(top, text=feature, font=LARGE_FONT).pack()
        Label(top,text="has no matches above the specified thresholds").pack()


def displayHit(result,feature,cluster,filename,ms1Deco,isoBounds):
    sample = next(result)
    resolution = 2
    top = Toplevel()
    Label(top,text=feature.replace("|","\n"),font=LARGE_FONT).pack()

    # w = 600
    # h = 800
    # x = 400
    # y = 100
    # top.geometry('%dx%d+%d+%d' % (w, h, x, y))
    #f = plt.figure(figsize = (20,10),dpi=100)
    f = plt.figure(figsize = (10,5),dpi=100)

    grid = plt.GridSpec(2,2,hspace=0.5,wspace=0.2)
    a1 = f.add_subplot(grid[0,:])
    a2 = f.add_subplot(grid[1,1])
    a3 = f.add_subplot(grid[1,0])

    indices = sample[3]
    decoSpec = sample[11]
    spec2Display = []
    ms1 = sample[12]

    for spec in sample[14:16]:
        tempSpec = np.zeros(MAXMASS*10**resolution)
        for i,m in zip(indices,spec):
            tempSpec[i] = m
        spec2Display.append(tempSpec.tolist())
    tempSpec = np.zeros(MAXMASS*10**resolution)
    for i,m in zip(indices,decoSpec):
        tempSpec[i] = m
    decoSpec = tempSpec.tolist()
    plotSpectra(0, maxMass=MAXMASS, res=0.01, spectras=spec2Display, specificPlot=a1,
                title="DP = " + str(np.round(sample[-7],2)) + "   PPM Error = " + str(np.round(((10**6)*(float(sample[2]) - float(sample[6]))/float(sample[2])),2)) + "    Redundant: " + str(sample[-1]),
                labels=["Component: " + sample[-4],sample[7]])
    spec = np.zeros(MAXMASS*10**resolution).tolist()
    for key,val in cluster["spectrum"].items():
        spec[int((10**resolution)*np.round(key,resolution))] += val
    plotSpectra(0, maxMass=MAXMASS, res=0.01, spectras=[spec,decoSpec], specificPlot=a2,
                title="DP = "+str(np.round(100*dotProductSpectra(spec,decoSpec),2)),
                labels=["Measured MS/MS Spectrum","Solved Spectrum"])
    ms1Spec = np.zeros(MAXMASS * 10 ** resolution).tolist()
    for m,i in ms1.items():
        ms1Spec[int((10 ** resolution) * np.round(m, resolution))] += i
    decoMs1Spec = np.zeros(MAXMASS * 10 ** resolution).tolist()
    for m,i in ms1Deco[sample[0]]:
        decoMs1Spec[int((10 ** resolution) * np.round(m, resolution))] += i

    a3.plot([isoBounds[0]-.5,isoBounds[1] + .5], [0, 0], c="black", linewidth=.5)
    plotSpectra(0, maxMass=MAXMASS, res=0.01, spectras=[ms1Spec, decoMs1Spec], specificPlot=a3,
                title="",
                labels=["Measured MS1 Spectrum", "Reconstructed Pure MS1 Spectrum"])
    canvas = FigureCanvasTkAgg(f, top)
    canvas._tkcanvas.pack()


    acceptMatchButton = ttk.Button(top, text="Accept Match",
                                 command=lambda: appendAnnotation(sample,cluster,filename))
    acceptMatchButton.pack()
    nextFigureButton = ttk.Button(top, text="Next",
                                 command=lambda: closeAndRestart(top,lambda:displayHit(result,feature,cluster,filename,ms1Deco,isoBounds)))
    nextFigureButton.pack()
    top.resizable = True


def appendAnnotation(sample,cluster,filename):
    if ".raw" in filename:
        filename = filename.replace(".raw","_Annotation.csv")
    elif ".mzML" in filename:
        filename = filename.replace(".mzML","_Annotation.csv")

    if not os.path.isfile(filename):
        file = open(filename,"w")
        file.write("#featureID,isolation_center_m/z,rt,compound_m/z,compound_rt,compound_formula,DB_Compound_ID,Compound_Name,DB_Spectrum_ID,dot_product,ppm_Error\n")
    else:
        file = open(filename,"a")

    toWrite = [cluster["group"],cluster["m/z"],sample[1],sample[6],sample[8],sample[9],sample[4],
               sample[7],sample[5],sample[13],1e6*abs(float(cluster["m/z"]) - float(sample[6]))/float(cluster["m/z"])]
    file.write(str(toWrite[0]))
    [file.write(","+str(x)) for x in toWrite[1:]]
    file.write("\n")
    file.close()

def writeResults(filenameOrig,scoreThresh,ppmThresh,usePPM,scanNumber,recursive,iso,peaks):

    try:
        if PPMFILTER[usePPM]:
            ppmThresh = np.inf
        if ".raw" in filenameOrig:
            ending = ".raw"
        else:
            ending = ".mzML"

        name = filenameOrig.replace(ending, "_decoID.csv")
        results = open(name,"r").readlines()
        filename = filenameOrig.replace(ending, "_Filtered.csv")
        file = open(filename, "w")
        file.write("#decoID Results for " + filenameOrig + " where the matches have a dot product > " + str(
            scoreThresh) + "and PPM <" + str(ppmThresh))
        file.write(results[0])

        for line in results[1:]:
            temp = line.rstrip().split(",")
            if scanNumber != "":
                if abs(float(temp[-3])) < ppmThresh and float(temp[-4]) > scoreThresh and temp[0] == scanNumber:
                    file.write(line)
            else:
                if abs(float(temp[-3])) < ppmThresh and float(temp[-4]) > scoreThresh:
                    file.write(line)
    except:
        top = Toplevel()
        Label(top, text="Datafile does not exist for selected file or is open in another program.\nPlease run the search first or close the csv file.",
              font=LARGE_FONT).pack()



class StartPage(Frame):


    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.Datafile = None
        #parent.config(bg="white")

        LeftTop = Frame(parent,width=500,height=500,pady=3,padx=3,bd=1)#,bg="white")

        label = Label(LeftTop, text="Select Parameters", font=LARGE_FONT)
        f = font.Font(label,label.cget("font"))
        f.configure(underline=True)
        label.configure(font=f)
        label.grid(row=0,column=0,pady=10, padx=10)

        perFragFrame = Frame(LeftTop)
        # Label(perFragFrame, text="% of Fragments Needed: ").grid(row=0,column=0,pady=10, padx=10)
        # perFrags = StringVar()
        # perFrags.set("0")
        # scaleFragment = Scale(perFragFrame, orient='horizontal', from_=0, to=100,command=lambda _:perFrags.set(scaleFragment.get()))
        # def val1(input):
        #     if input == "":
        #         #scale.set(0)
        #         return True
        #     try:
        #         t = int(input)
        #         if t >= 0 and t <= 100:
        #             scaleFragment.set(int(input))
        #             return True
        #         else:
        #             return False
        #     except:
        #         return False

        # e1 = Entry(perFragFrame, width=3, textvariable=perFrags, validate="key",
        #           vcmd=(perFragFrame.register(val1), ('%P',)))
        # e1.grid(row=0, column=1, pady=10, padx=10)
        #
        # scaleFragment.grid(row=0,column=2,pady=10, padx=10)
        # perFragFrame.grid(row=1,column = 0)


        procNumFrame = Frame(LeftTop)
        Label(procNumFrame, text="Processor Number").grid(row=0,column=0,pady=10, padx=10)
        numCores = StringVar()
        numCores.set("2")

        menu = OptionMenu(procNumFrame, numCores, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
        menu.grid(row=0,column=1,pady=10, padx=10)
        procNumFrame.grid(row=1,column=0)

        useRT = Radiobar(LeftTop, ['Yes', 'No'], label="Use Retention Time:")
        useRT.grid(row=2,column=0,pady=10, padx=10)

        RTFrame = Frame(LeftTop)
        Label(RTFrame,text="Retention Time Tolerance (minutes)").grid(row=0,column=0,pady=10)
        rtTol = StringVar()
        rtTol.set(".5")
        Label(RTFrame,textvariable = rtTol).grid(row=0,column=1,pady=10,padx=10)

        def val3(input):
            if input == "": return True
            try:
                t = float(input)
                rtTol.set(str(t))
                return True
            except:
                isoWindow.set("")
                return False

        e = Entry(RTFrame, width=3, textvariable=rtTol, validate="key",
                  vcmd=(RTFrame.register(val3), ('%P',)))
        e.grid(row=0, column=2, pady=10, padx=10)

        RTFrame.grid(row=3,column=0,pady=10,padx=10)

        standardOrDeco = Radiobar(LeftTop, ['Pairwise', 'Deconvolution'], label="Search Method:")
        standardOrDeco.grid(row=4,column=0,pady=10, padx=10)

        recursive = Radiobar(LeftTop, ['Yes', 'No'], label="Use Predicted Unknown Library:")
        recursive.grid(row=5,column=0,pady=10, padx=10)

        iso = Radiobar(LeftTop, ['Yes', 'No'], label="Use Simulated M+1 Isotopologue Spectra:")
        iso.grid(row=6,column=0,pady=10, padx=10)

        peaks = Radiobar(LeftTop, ['Yes', 'No'], label="Use MS1 data:")
        peaks.grid(row=7,column=0,pady=10, padx=10)

        dtype = Radiobar(LeftTop, ['DDA', 'DIA'], label="Data Acquisition Method:")
        dtype.grid(row=8,column=0,pady=10, padx=10)

        massAccFrame = Frame(LeftTop)
        Label(massAccFrame, text="Mass PPM Tolerance").grid(row=0,column=0,pady=10, padx=10)
        scaleAcc = Scale(massAccFrame, orient='horizontal', from_=1, to=50)
        scaleAcc.set(5)
        scaleAcc.grid(row=0,column=1,pady=10, padx=10)
        massAccFrame.grid(row=9,column=0)

        fragIntFrame = Frame(LeftTop)
        Label(fragIntFrame,text="Fragment Intensity Threshold").grid(row=0,column=0,pady=10)
        fragInt = StringVar()
        fragInt.set("0")
        Label(fragIntFrame,textvariable = fragInt).grid(row=0,column=1,pady=10,padx=10)

        def val3(input):
            if input == "": return True
            try:
                t = float(input)
                fragInt.set(str(t))
                return True
            except:
                fragInt.set("")
                return False

        e = Entry(fragIntFrame, width=3, textvariable=fragInt, validate="key",
                  vcmd=(fragIntFrame.register(val3), ('%P',)))
        e.grid(row=0, column=2, pady=10, padx=10)

        fragIntFrame.grid(row=10,column=0)

        Label(LeftTop,text="Isolation Window Width (Da, required for non-thermo data)").grid(row=11,column=0)
        isoWindowFrame = Frame(LeftTop)
        isoWindow = StringVar()
        isoWindow.set("1")
        isoWindow.set("1")

        Label(isoWindowFrame,textvariable = isoWindow).grid(row=0,column=1,pady=10,padx=10)

        def val3(input):
            if input == "": return True
            try:
                t = float(input)
                isoWindow.set(str(t))
                return True
            except:
                isoWindow.set("")
                return False
        e = Entry(isoWindowFrame, width=3, textvariable=isoWindow, validate="key",
                  vcmd=(isoWindowFrame.register(val3), ('%P',)))
        e.grid(row=0, column=2, pady=10, padx=10)
        isoWindowFrame.grid(row=12, column=0)


        RightTop = LeftTop#Frame(parent,width=1000,height=1000,pady=3,padx=3,bd=1)#,bg="white")

        lab = Label(RightTop, text="Select MS/MS Library", font=LARGE_FONT)
        lab.grid(row=0,column = 3,pady=10,padx=10)
        f = font.Font(label,label.cget("font"))
        f.configure(underline=True)
        lab.configure(font=f)

        custOrMzCloud = Radiobar(LeftTop, ['mzCloud (reference)','mzCloud (autoprocessing)', "MoNA", "HMDB",'Custom Library'], label="Library Choice:")
        custOrMzCloud.grid(row=1,column=3)

        libraryFile = StringVar()
        libraryFile.set("No File Selected")
        Label(LeftTop,textvariable=libraryFile).grid(row=2,column=3,pady=10)
        browseLibFileButton = ttk.Button(LeftTop, text="Select File",
                                      command=lambda: browseForFile(libraryFile,types=[("database",".db"),("msp",".msp")]))
        browseLibFileButton.grid(row=3,column=3,pady=10)


        mzCloudAPIFile = StringVar()
        mzCloudAPIFile.set("No File Selected")
        Label(LeftTop,textvariable=mzCloudAPIFile).grid(row=4,column=3,pady=10)
        browsemzCloudKeyFileButton = ttk.Button(LeftTop,text="Select mzCloud API Key File",
                                                command=lambda:browseForFile(mzCloudAPIFile,types=[("txt","*.txt")]))
        browsemzCloudKeyFileButton.grid(row=5,column=3,pady=10)

        lab = Label(RightTop, text="Select MS Datafile", font=LARGE_FONT)
        lab.grid(row=6,column = 3,pady=10,padx=10)
        f = font.Font(label,label.cget("font"))
        f.configure(underline=True)
        lab.configure(font=f)

        filename = StringVar()
        filename.set("No File Selected")
        Label(LeftTop, textvariable=filename).grid(row=7,column=3,pady=10)
        browseFileButton = ttk.Button(LeftTop, text="Select File",
                                      command=lambda: browseForFile(filename,types=[("mzML",".mzML"),("Thermo", "*.raw"),("Agilent","*.d"),("DecoID (visualization only)","*.DecoID"),("All Files","*")]))
        browseFileButton.grid(row=8,column=3,pady=10)

        lab = Label(RightTop, text="Select Peak Information File", font=LARGE_FONT)
        lab.grid(row=9, column=3, pady=10, padx=10)
        f = font.Font(label, label.cget("font"))
        f.configure(underline=True)
        lab.configure(font=f)

        filenamePeak = StringVar()
        filenamePeak.set("No File Selected")
        Label(LeftTop, textvariable=filenamePeak).grid(row=10, column=3, pady=10)
        browsePeakFileButton = ttk.Button(LeftTop, text="Select File (optional for DDA)",
                                      command=lambda: browseForFile(filenamePeak,
                                                                    types=[("csv", "*.csv"),
                                                                           ("All Files", "*")]))
        browsePeakFileButton.grid(row=11, column=3, pady=10)
        runSearchButton = ttk.Button(LeftTop, text="Search",
                                     command=lambda: performSearch(filename,numCores,recursive,iso,peaks,dtype,scaleAcc,custOrMzCloud,libraryFile,standardOrDeco,filenamePeak,float(isoWindow.get())/2,useRT,rtTol,mzCloudAPIFile,fragInt))
        runSearchButton.grid(row=12,column = 3,pady
        =10, padx=10)


        LeftTop.grid(row=0,column=0)




        lab = Label(RightTop, text="Analyze Search Results", font=LARGE_FONT)
        lab.configure(font=f)
        lab.grid(row=0,column=6,pady=10, padx=10)

        scoreThreshFrame = Frame(RightTop)
        Label(scoreThreshFrame, text="Score Threshold: ").grid(row=0,column=0,pady=10, padx=10)
        scoreThresh = StringVar()
        scoreThresh.set("50")
        scale = Scale(scoreThreshFrame, orient='horizontal', from_=0, to=100,command=lambda _:scoreThresh.set(scale.get()))
        scale.set(50)
        def val(input):
            if input == "":
                #scale.set(0)
                return True
            try:
                t = int(input)
                if t >= 0 and t <= 100:
                    scale.set(int(input))
                    return True
                else:
                    return False
            except:
                return False
        e = Entry(scoreThreshFrame,width=3,textvariable=scoreThresh,validate="key",vcmd=(scoreThreshFrame.register(val),('%P',)))
        e.grid(row=0,column=1,pady=10, padx=10)
        scale.grid(row=0,column=2)
        scoreThreshFrame.grid(row=1,column=6)


        PPMThreshold = Radiobar(RightTop, ['Yes', 'No'], label="Use PPM Threshold:")
        PPMThreshold.grid(row=2,column=6,pady=10, padx=10)

        ppmThreshFrame = Frame(RightTop)
        Label(ppmThreshFrame, text="PPM Threshold").grid(row=0,column=0,pady=10, padx=10)
        scalePPM = Scale(ppmThreshFrame, orient='horizontal', from_=1, to=50)
        scalePPM.grid(row=0,column=1,pady=10, padx=10)
        ppmThreshFrame.grid(row=3,column=6)

        scanNumberFrame = Frame(RightTop)
        Label(scanNumberFrame, text="Feature ID: ").grid(row=0,column=0,pady=10, padx=10)
        scanNumber = StringVar()
        scanNumber2 = StringVar()
        scanNumber.set("")
        scanNumber2.set("")

        Label(scanNumberFrame,textvariable = scanNumber2).grid(row=0,column=1,pady=10,padx=10)

        def val2(input):
            if input == "": return True
            try:
                t = int(input)
                scanNumber2.set(str(t))
                return True
            except:
                scanNumber2.set("")
                return False
        e = Entry(scanNumberFrame, width=6, textvariable=scanNumber, validate="key",
                  vcmd=(scanNumberFrame.register(val2), ('%P',)))
        e.grid(row=0, column=2, pady=10, padx=10)
        scanNumberFrame.grid(row=4, column=6)


        visualizeResultsButton = ttk.Button(RightTop, text="Display Search Results",
                                            command=lambda: visualizeResultsStart(self,filename,scale.get(),scalePPM.get(),PPMThreshold.state(),scanNumber2.get(),recursive,iso,peaks,dtype.state()))
        visualizeResultsButton.grid(row=5,column=6,pady=10, padx=10)

        writeResultsButton = ttk.Button(RightTop, text="Write Filtered Results",
                                        command=lambda: writeResults(filename.get(),scale.get(),scalePPM.get(),PPMThreshold.state(),scanNumber2.get(),recursive,iso,peaks))
        writeResultsButton.grid(row=6,column=6,pady=10, padx=10)
