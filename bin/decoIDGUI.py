import sys
from ctypes import windll
import os
import multiprocessing

if __name__ == '__main__':
    multiprocessing.freeze_support()
    multiprocessing.set_start_method("spawn")
    #windll.shcore.SetProcessDpiAwareness(1)
    sys.path.append(os.path.join(os.path.dirname(__file__),"..","src/"))

    from GUIHelper import *

    app = decoIDSearch()
    app.mainloop()