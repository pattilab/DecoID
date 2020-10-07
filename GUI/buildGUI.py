import PyInstaller.__main__
import os



PyInstaller.__main__.run([
    '--name=%s' % "DecoIDGUI",
    '--onedir',
    "--hidden-import=sklearn.neighbors.typedefs",
    "--hidden-import=sklearn.neighbors.quad_tree",
    "--hidden-import=sklearn.tree._utils",    #'--windowed',
    "--hidden-import=sklearn.utils._cython_blas",
    "--hidden-import=DecoID.DecoID",
    "--additional-hooks-dir=hook-sklearn.linear_model.py",
    '--add-binary=%s' % 'C:/Users/Ethan/AppData/Local/Programs/Python/Python37/Lib/site-packages/sklearn/.libs/vcomp140.dll;.', #need to change to local directory of the .dll file
    '--add-binary=%s' % "../databases/HMDB_experimental.db;.",
    '--add-binary=%s' % "../databases/MoNA-export-LC-MS-MS_Spectra.db;.",
    '--add-binary=%s' % '../src/DecoID/mzCloudCompound2TreeLinkageautoprocessing.pkl;.',
    '--add-binary=%s' % '../src/DecoID/mzCloudCompound2TreeLinkagereference.pkl;.',

    '-y',
    'decoIDGUI.py'])

