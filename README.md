# DecoID
Metabolomics software for database-assisted deconvolution of MS/MS spectra

## System Requirements

Standalone executible was built for Windows 10 64 bit

Package has been tested with Python 3.7

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
libFile = "DecoID/databases/HMDB_experimental.db" #path to database
numCores = 10 # # of parallel processes to use
file = "DecoID/exampleData/Asp-Mal_1uM_5Da.mzML" #path to datafile
peakfile = "DecoID/exampleData/peak_table.csv" #path to peak information file

useMS1 = True #use MS1 data if avaible
DDA = True #Data is DDA
massAcc = 10 #Mass accuracy of instrument
useIso = True #remove contamination from orphan isotopologues
lam = 1 #LASSO regularizing coefficient (recommendation: 1 for DDA, 100 for DIA)

if __name__ == '__main__':
    decID = DecoID(libFile,numCores) #create DecoID object
    decID.readData(file, 2, useMS1, DDA, massAcc,peakDefinitions=peakfile) #load in datafile
    decID.identifyUnknowns() #identify compounds for inclusion in on-the-fly unknown library
    decID.searchSpectra("y", lam , iso = useIso) #deconvolve and identify spectra

```

Example usage on DIA MS/MS datafile (larger and slower, >20 min)

```
from DecoID.DecoID import DecoID
libFile = "DecoID/databases/HMDB_experimental.db" #path to database
numCores = 10 # # of parallel processes to use
file = "DecoID/exampleData/IROA_P1-6_DIA_test_pos1.mzML" #path to datafile
peakfile = "DecoID/exampleData/IROA_p1-6_peak_table_pos_v3.csv" #path to peak information file

useMS1 = True #use MS1 data if avaible
DDA = False #Data is DIA
massAcc = 10 #Mass accuracy of instrument
useIso = False # do not remove contamination from orphan isotopologues (assume negligent for DIA)
lam = 100 #LASSO regularizing coefficient (recommendation: 1 for DDA, 100 for DIA)

if __name__ == '__main__':
    decID = DecoID(libFile,numCores) #create DecoID object
    decID.readData(file, 2, useMS1, DDA, massAcc,peakDefinitions=peakfile) #load in datafile
    decID.searchSpectra("y", lam , iso = useIso) #deconvolve and identify spectra
```

expected output files are included in the exampleData directory




