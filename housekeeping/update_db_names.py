import pickle as pkl

link = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkage_InChIreference.pkl", "rb"))

linkerKey = {pol:{x[2]:x for x in link[pol]} for pol in link}
db = pkl.load(open("../databases/mzCloud_reference.db","rb"))
newDB = {}

for pol in db:
    print(pol)
    newDB[pol] = {}
    for id,spec in db[pol].items():
        name = spec["name"]
        form = spec["formula"]
        if "M+1" in name:
            m1 = True
            name = name.replace(" (M+1)","")
            addition = " (M+1)"
        else:
            m1 = False
            addition = ""
        info = linkerKey[pol][name]
        newDB[pol][id] = spec
        newDB[pol][id]["cpdID"] = info[1] + addition

        if len(newDB[pol]) % 100:
            print(len(newDB[pol])/len(db[pol]))

pkl.dump(newDB,open("../databases/mzCloud_reference_InChIKey.db","wb"))

