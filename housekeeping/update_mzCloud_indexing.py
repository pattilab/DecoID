from DecoID.DecoID import DecoID,mzCloudPy
import pickle as pkl

api_key = open("mzCloud_api_key.txt","r").readline().rstrip()

#mzCloudPy.generateCompoundID2SpectralIDIndexedByM_ZStrict(100,api_key,"reference")
#mzCloudPy.generateCompoundID2SpectralIDIndexedByM_ZStrict(100,api_key,"autoprocessing")

tmp = pkl.load(open("../src/DecoID/mzCloudCompound2TreeLinkage_InChI" + "reference" + ".pkl", "rb"))

print(len(tmp["Positive"]))

for x in tmp["Positive"]:
    print(x)
    break