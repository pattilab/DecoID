from DecoID.DecoID import DecoID,mzCloudPy

api_key = open("mzCloud_api_key.txt","r").readline().rstrip()

#mzCloudPy.generateCompoundID2SpectralIDIndexedByM_ZStrict(100,api_key,"reference")
mzCloudPy.generateCompoundID2SpectralIDIndexedByM_ZStrict(100,api_key,"autoprocessing")