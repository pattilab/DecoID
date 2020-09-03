from DecoID.DecoID import DecoID,mzCloudPy,Keys
import pickle as pkl

library = "autoprocessing"

if library == "reference":
    link = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkagereference.pkl","rb"))
    prefix = "r"
else:
    link = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkageautoprocessing.pkl","rb"))
    prefix = "a"

api_key = open("mzCloud_api_key.txt","r").readline().rstrip()
key = Keys(api_key)
db = {pol:{} for pol in link}
for pol in link:
    print(len(link[pol]))
    ind = 0
    for tree in link[pol]:
        specTree = mzCloudPy.getTrees([tree],key)
        ind += 1
        for specID in specTree:
            db[pol][specID] = {"cpdID":prefix+str(tree[1]),"id":specID,"rt":tree[4],"formula":tree[3],"name":tree[2].replace(",","_"),"mode":pol,"spec":specTree[specID],"m/z":tree[5]}
        if ind % 100 == 0:
            print(ind/len(link[pol]))
pkl.dump(db,open("../databases/mzCloud_"+library+".db","wb"))
