from DecoID.DecoID import DecoID,mzCloudPy,Keys,getSubForms,createM1FromM0andFragAnnotation
import pickle as pkl
from multiprocessing import Queue,Process
from copy import deepcopy

# library = "autoprocessing"
#
# if library == "reference":
#     link = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkagereference.pkl","rb"))
#     prefix = "r"
# else:
#     link = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkageautoprocessing.pkl","rb"))
#     prefix = "a"
mPlusOnePPM = 15
numCores = 40
# api_key = open("mzCloud_api_key.txt","r").readline().rstrip()
# key = Keys(api_key)
# db = {pol:{} for pol in link}
# for pol in link:
#     print(len(link[pol]))
#     ind = 0
#     for tree in link[pol]:
#         specTree = mzCloudPy.getTrees([tree],key)
#         ind += 1
#
#         if len(specTree) > 0:
#             for specID in specTree[tree]:
#                 db[pol][specID] = {"cpdID":prefix+str(tree[1]),"id":specID,"rt":tree[4],"formula":tree[3],"name":tree[2].replace(",","_"),"mode":pol,"spec":specTree[tree][specID],"m/z":tree[5]}
#         if ind % 100 == 0:
#             print(ind/len(link[pol]))

def processCompoundGroup(cpd,spectraForCompounds,pol,q):
    specIDs = list(spectraForCompounds.keys())
    formula = spectraForCompounds[specIDs[0]]["formula"]
    cpdID = str(cpd) + " (M+1)"
    name = spectraForCompounds[specIDs[0]]["name"] + " (M+1)"
    masses, sub_forms, carbon_pos, bounds = getSubForms(formula, pol)
    for specID in spectraForCompounds:
        key = str(specID) + "_M1"
        val = deepcopy(spectraForCompounds[specID])
        val["id"] = key
        val["cpdID"] = cpdID
        val["name"] = name
        val["m/z"] = val["m/z"] + 1.00335
        val["spec"],_ = createM1FromM0andFragAnnotation(val["spec"],masses,sub_forms,carbon_pos,bounds,mPlusOnePPM)
        q.put([key,val,pol])

library = "reference"
if __name__ == '__main__':

    db = pkl.load(open("../databases/mzCloud_reference.db","rb"))
    db = {pol:{key:val for key,val in db[pol].items() if "M+1" not in val["cpdID"]} for pol in db}
    q = Queue()
    processes = []
    cpdsDict = {}
    for pol in db:
        cpdsDict[pol] = {}
        for specID in db[pol]:
            if db[pol][specID]["cpdID"] not in cpdsDict[pol]:
                cpdsDict[pol][db[pol][specID]["cpdID"]] = {}
            cpdsDict[pol][db[pol][specID]["cpdID"]][specID] = dict(db[pol][specID])
    started = 0
    for pol in db:
        print(pol)
        for cpd in cpdsDict[pol]:
            p = Process(target=processCompoundGroup, args=(cpd,cpdsDict[pol][cpd],pol,q))
            while len(processes) >= numCores:
                if not q.empty():
                    k,v,m = q.get()
                    if len(v["spec"]) > 0:
                        db[m][k] = v
                processes = [x for x in processes if x.is_alive()]
            p.start()
            started += 1
            if started % 100 == 0:
                print(started/len(cpdsDict[pol]))
            processes.append(p)

    while len(processes) > 0:
        if not q.empty():
            k, v, m = q.get()
            if len(v["spec"]) > 0:
                db[m][k] = v
        processes = [x for x in processes if x.is_alive()]
    while not q.empty():
        k, v, m = q.get()
        if len(v["spec"]) > 0:
            db[m][k] = v

    pkl.dump(db,open("../databases/mzCloud_"+library+"_correct.db","wb"))
