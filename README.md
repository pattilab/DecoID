# DecoID
Metabolomics software for database-assisted deconvolution of MS/MS spectra

## System Requirements

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

API Documentation: 

## Demo

Demo data available under DecoID/exampleData/

Example Usage of Deconvolution of DDA datafile:

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
    decID.readData(file, 2, usePeaks, DDA, massAcc,peakDefinitions=peakfile) #load in datafile
    decID.identifyUnknowns() #identify compounds for inclusion in on-the-fly unknown library
    decID.searchSpectra("y", lam , iso = useIso) #deconvolve and identify spectra
```
