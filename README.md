# DecoID
Metabolomics software for database-assisted deconvolution of MS/MS spectra

## System Requirements

Standalone executible was built for Windows 10 64 bit

Package has been tested with Python 3.7 on Windows 10 and macOS Catalina 

Package has the following dependencies:

numpy (v1.18.1)

sklearn (v0.22.1)

pandas (v1.0.1)

dill (v0.3.1.1)

scipy (v1.4.1)

pyteomics (v4.2)

requests (v2.22.0)

lxml (v4.5.0)

molmass (2020.6.10)

keras (v2.4.3)

tensorflow (v2.4.0)

IsoSpecPy (v2.1.4)

Memory usage can be very intensive when searching DIA data or MS/MS spectra acquired with wide isolation windows (>10 m/z). This may limit the number of parallel processes
that can be run. If memory errors occur, reduce this value. 


In order to process vendor formatted data without manual conversion, MS-Convert (http://proteowizard.sourceforge.net/tools.shtml) needs to be installed and added to PATH. 

## Installation

### Installation with ```pip```:

```
pip install DecoID
```
PyPI:

https://pypi.org/project/DecoID/

### Manual installation from source:

```
git clone https://github.com/e-stan/DecoID.git
pip install DecoID/src/.
```

### Installation of Standalone Windows .exe for User Interface

Download zip from http://pattilab.wustl.edu/software/DecoID

unzip file and run DecoIDGUI/DecoIDGUI.exe

## Usage

User interface documentation and guide available at DecoID/DecoIDGUI_manual.pdf

Total size is approximatly 1gb and includes binaries of HMDB and MoNA. Installation time is dependent on network speed. With around 100 Mb/sec download speed, total time to download/extract/run was approximatly 5 minutes.

API Documentation: https://decoid.readthedocs.io/

## Demo

Demo data available under DecoID/exampleData/

Example usage of deconvolution of DDA datafile (fast <5min):

```
from DecoID.DecoID import DecoID

#sets database to use
libFile = "../databases/HMDB_experimental.db"

#mzCloud key if necessary
key = "none"
mzCloudLib = "reference"

#number of parallel processes to use
numCores = 4

#filename of query MS/MS data
file = "../exampleData/Asp-Mal_1uM_5Da.mzML"

#filename of peak list
peakFile = "../exampleData/peak_table.csv"

#set parameters
usePeaks = True
DDA = True #data is DDA
massAcc = 10 #ppm tolerance
fragThresh= 0.01 #require non-zero dot product threshold
offset = .5 #half of isolation window width. Only for non-thermo data
useIso = True #use predicted M+1 isotopolgoue spectra
threshold = 0 #minimum dot product for reporting
lam = 5.0 #LASSO parameter
rtTol = float("inf") #retention time tolerance for database, inf means ignore RT
fragCutoff = 1000 #intensity threshold for MS/MS peaks


if __name__ == '__main__':

    #create DecoID object
    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key)

    #read in data
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions=peakFile,frag_cutoff=fragCutoff)

    #identify unknowns compounds for on-the-fly unknown library
    decID.identifyUnknowns(iso=useIso,rtTol=rtTol,dpThresh=80,resPenalty=lam)

    #search spectra
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)

```

Example usage on DIA MS/MS datafile (larger and slower, >20 min)

```
from DecoID.DecoID import DecoID

#sets database to use
libFile = "../databases/HMDB_experimental.db"

#mzCloud key if necessary
key = "none"
mzCloudLib = "reference"

#number of parallel processes to use
numCores = 5

#filename of query MS/MS data
file = "../exampleData/IROA_P1-6_DIA_test_pos1.mzML"

#filename of peak list
peakFile = "../exampleData/IROA_p1-6_peak_table_pos_v3.csv"

#set parameters
usePeaks = True
DDA = False #data is DIA
massAcc = 10 #ppm tolerance
fragThresh= 0.01 #require non-zero dot product threshold
offset = .5 #half of isolation window width. Only for non-thermo data
useIso = True #use predicted M+1 isotopolgoue spectra
threshold = 0 #minimum dot product for reporting
lam = 50.0 #LASSO parameter
rtTol = float("inf") #retention time tolerance for database, inf means ignore RT
fragCutoff = 1000 #intensity threshold for MS/MS peaks

if __name__ == '__main__':

    #create DecoID object
    decID = DecoID(libFile, mzCloudLib, numCores,api_key=key)

    #read in data
    decID.readData(file, 2, usePeaks, DDA, massAcc,offset,peakDefinitions=peakFile,frag_cutoff=fragCutoff)

    #search spectra
    decID.searchSpectra("y", lam , fragThresh, useIso, threshold,rtTol=rtTol)

```

expected output files are included in the exampleData directory




