DecoID Usage
==================================

Input Files
-----------

MS/MS datafile
++++++++++++++
DecoID supports all vendor file formats
compatible with MS-Convert. Natively, DecoID accepts
.mzML files that have been centroided. However, if MS-Convert
is installed on the machine and msconvert.exe is in the system
PATH, vendor files can be supplied. DecoID will read MS/MS metadata including polarity, targeted m/z's
, and isolation windows (Thermo Fisher Scientific datafiles only).
For non-Thermo data, the width of the isolation window must be
provided. DecoID is compatable with multiple MS/MS datatypes including
data dependent acquistion (DDA) and data independent acqusition (DIA). Further,
less common experimental workflows where MS1 full scan data is not collected can
also be used.

Peak Information
++++++++++++++++
To facilitate metabolite information and improve the signal to noise of the MS/MS
spectra and the resulting deconvolution, the m/z and retention time bounds of features
of interest can be provided in a .csv file. If so, DecoID will average all acquired spectra within the
retention time bounds. In addition, the simplifies the output datafile for downstream
analysis. Peak information is required for DIA data.

The format of this file is shown below:

======   ========   =====
mz       rt_start   rt_end
======   ========   =====
123.45   1.67       2.43
145.78   2.94       3.55
90.08    7.83       9.23
etc.     etc.       etc.
======   ========   =====

Database
++++++++
DecoID uses either NIST msp formatted databases or a remote connection to the mzCloud
(Thermo-Fisher Scientific) spectral database to perform metabolite identification and
MS/MS deconvolution. For mzCloud based deconvolution, no input file is required. More
details are given below in `Connection to mzCloud`_. MSP files of the Human Metabolome Database (HMDB)
and the Mass Bank of North America and be downloaded from the `DecoID GitHub Site <https://github.com/e-stan/DecoID/releases>`_
Alternatively, these databases can be downloaded from the MoNA website `here <https://mona.fiehnlab.ucdavis.edu/downloads>`_.
MSP files are simple text files. For compatibility with DecoID the following fields must be non-empty:

* Name
* DB#
* Precursor_type OR Ion_mode
* ExactMass OR PrecursorMZ

An example MSP entry is shown below:::

    Name: GLYCEROL
    Synon: $:00in-source
    DB#: 0
    InChiKey: 56-81-5
    Precursor_type: [M-H]-
    Spectrum_type: MS2
    PrecursorMZ: 91.04004
    Ion_mode: N
    Collision_energy: HCD (NCE 40%)
    RetentionTime: 2.969
    Formula: C3H8O3
    MW: 92.0
    ExactMass: 92.04734
    Num Peaks: 40
    50.003 0.002940881981667288
    50.043 0.001125309959124438
    52.92 0.0010507693668541853
    54.262 0.0009862185551523176
    54.741 0.0010954343864214948
    55.606 0.0012054987503008479
    55.851 0.0010793672216122406
    58.039 0.0010972036230296137
    58.056 0.0012052440367981193
    58.173 0.0011118376593017905
    59.013 0.1877307037191327
    59.017 0.004254219849821239
    59.967 0.9227487450965145
    59.97 0.038360551398232204
    59.972 0.023916371530628082
    60.974 0.0013073255841957638
    61.744 0.001150458036919167
    62.524 0.0010455834798728655
    62.724 0.001173758557930979
    63.886 0.0012296895388148035
    64.873 0.0010963706929665274
    65.145 0.001130812508637816
    72.444 0.0012381777855903992
    72.449 0.001096971813758573
    73.011 0.0046482667112914475
    74.99 1.0
    74.995 0.03924297796715474
    74.997 0.026938366619882102
    75.967 0.0011636634025248632
    76.214 0.0012321326056216008
    76.969 0.005576439025961694
    91.021 0.002564739158271785
    91.029 0.0013610635233254089
    93.001 0.9799805233973448
    93.007 0.027811413214189975
    93.01 0.004905488457969213
    93.011 0.004496375838517466
    97.159 0.0011219222849101172
    98.845 0.0012869953884785937
    100.86 0.0010522699014690819

After DecoID parses the MSP file a binary .db file is generated (this is a Pickle file) for faster loading
on future usage.

Output Files
------------
After successful deconvolution of an MS/MS datafile, 3 output files are generated in the input file directory.

* <fn>_scanInfo.csv
* <fn>_decoID.csv
* <fn>.DecoID

**<fn>_scanInfo.csv** gives the purified spectra and all components for each acquired MS/MS spectrum. It is formatted
as shown below:

======  =====================   =============   =========== ================== ============ ==== ==================================
scanID  Signal to Noise Ratio   numComponents   componentID componentAbundance  ComponentMz rt   spectrum
======  =====================   =============   =========== ================== ============ ==== ==================================
1       2.5                     4               cpdID1      .8                 133.014      2.4  84:12.6 72:45
1       2.5                     4               cpdID2      .2                 132.54       2.4  102:45 68:55
1       2.5                     4               original    0                  133.014      2.4  68:55 72:45 84:12.6 102:45 72:55
1       2.5                     4               residual    0                  133.014      2.4  72:10
etc     etc                     etc             etc         etc                etc          etc  etc
======  =====================   =============   =========== ================== ============ ==== ==================================

* scanID: The row number from the peak information file that gives which feature this spectrum belongs or if peak information is not provided, this is the scanID of the MS/MS spectrum.
* Signal to Noise Ratio: denotes the  the signal attributed to a particular compound divided by the signal not attributed to any compound.
* numComponents: The number of components used in the deconvolution.
* componentID: compound ID for each component. If this field is "original" then this is the acquired spectrum. Residual is the error between the reconstructed spectrum and the acquired.
* ComponentMz: The m/z value of the precursor of the component.
* rt: The retention time where the spectrum was acquired.
* spectrum: The spectrum for each componet given in a m/z:intensity pairs separated by a space.


**<fn>_decoID.csv** gives the metabolite identification results after the deconvolution. With a single match on each line. The format is given below:

======  ====================    ====  ============   ===============  =============   ==============   ===========    =========   =========   ===========
scanID  isolation_center_m/z    rt    compound_m/z    DB_Compound_ID  Compound_Name   DB_Spectrum_ID   dot_product    ppm_Error   Abundance   ComponentID
======  ====================    ====  ============   ===============  =============   ==============   ===========    =========   =========   ===========
1       133.014                 2.4   133.014        cpdID01          Malic Acid      HMDB0031518      99.8           -1.3        0.8         cpdID01
etc     etc                     etc   etc            etc              etc             etc              etc            etc         etc         etc
======  ====================    ====  ============   ===============  =============   ==============   ===========    =========   =========   ===========

* scanID: The row number from the peak information file that gives which feature this spectrum belongs or if peak information is not provided, this is the scanID of the MS/MS spectrum.
* isolation_center_m/z: The feature of interest m/z value.
* rt: The retention time where the spectrum was acquired.
* compound_m/z: The m/z value of the matched compound.
* DB_Compound_ID: The compound ID of the matched compound.
* Compound_Name: The name of the matched compound.
* DB_Spectrum_ID: Spectrum ID or accession of the matched spectrum. Given by DB# in the input database.
* dot_product: The normalized dot product similarity to the reference spectrum
* ppm_Error: The mass error in parts per million (ppm) between the feature's m/z and the database match m/z.
* Abundance: The normalized regression coefficient of this compound in the deconvolution. Note: this should not be used for comparative/quantitative purposes.
* componentID: compound ID of the component matched to. If this field is "original" then this is the acquired spectrum. Residual is the error between the reconstructed spectrum and the acquired.

**<fn>.DecoID** is a gzipped pickle file that contains all the information provided in the previous two output files but in a format that allows for easier analysis and visualation through
the DecoID user interface.

Example Usage
--------------------

Regardless of data type the following parameters are required:::

    from DecoID.DecoID import DecoID
    libFile = "DecoID/databases/HMDB_experimental.db" #path to database
    numCores = 10 # of parallel processes to use
    file = "DecoID/exampleData/Asp-Mal_1uM_5Da.mzML" #path to datafile
    peakfile = "DecoID/exampleData/peak_table.csv" #path to peak information file

    useMS1 = True #use MS1 data if available
    massAcc = 10 #Mass accuracy of instrument
    res = 2 # # of decimal places to round MS/MS peaks.

With this the DecoID object can be instantiated and database parsed:::

    decID = DecoID(libFile,numCores)

Before the raw data can be read-in some data-type specific parameters must be provided:::

    offset = .5 #half of the width of the MS/MS isolation window. Not required for Thermo data.
    DDA = True #true for DDA, False for DIA

Now the raw MS/MS data can be read:::

    decID.readData(file, 2, useMS1, DDA, massAcc,peakDefinitions=peakfile)

With the data read, the search parameters can be defined:::

    lam = 1 # LASSO regression coefficient. The higher this is the more sparse a solution will be found. Recommend 1 for DDA and 100 for DIA.
    useIso = True # Predict M+1 isotopologue spectra to remove contamination from orphan isotopologues

Optionally, acquired pure MS/MS spectra can be used to deconvolve spectra in the datafile if data is from a DDA experiment. To enable this the command below must be run:::

    decID.identifyUnknowns(resPenalty=lam,iso=useIso)

Now, the datafile can be searched with the command below:::

    decID.searchSpectra("y",lam,iso=useIso)

Advanced Usage
--------------

Changing Deconvolution Parameters
+++++++++++++++++++++++++++++++++

Changing the "lam" parameter in effect allows for a continuum of performance between direct library searching without deconvolution and non-regularized deconvolution.

With::

    lam = float("inf")

You have no deconvolution and standard library searching.

With::

    lam = 0

There is no penalty for more complex solutions.

High Performance Computing
++++++++++++++++++++++++++

Documentation in progress. See DecoID/HPC_scripts/

Connection to mzCloud
+++++++++++++++++++++

Connection to mzCloud is dependent on an access key granted by Thermo-Fisher Scientific. If a key is granted, it must
be entered during instantiation of the DecoID object and the libFile parameter must be "none":::

    decID = DecoID("none",numCores,api_key="XXXXXXXXXXXX")


