from scipy.interpolate import griddata
import numpy as np
import time
from matplotlib import pyplot

import concurrent.futures
import time



class Test:
    def __init__(self,length):
        self.list = list(np.zeros(length))

    def do_something(self, seconds):
        print(self.list[seconds])

start = time.perf_counter()
a = Test(5)
a.list[1] = 2
a.do_something(2)




if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        secs = [5, 4, 3, 2, 1]
        results = executor.map(a.do_something, secs)

        # for result in results:
        #     print(result)

    finish = time.perf_counter()

    print(f'Finished in {round(finish-start, 2)} second(s)')


