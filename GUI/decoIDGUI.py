
import multiprocessing

if __name__ == '__main__':
    multiprocessing.freeze_support()
    multiprocessing.set_start_method("spawn")

    import GUIHelper

    from GUIHelper import *

    app = decoIDSearch()
    app.mainloop()