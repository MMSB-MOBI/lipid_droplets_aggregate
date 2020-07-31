import time
import numpy as np
import multiprocessing as mp

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[2000000, 5])
data = arr.tolist()

start = time.time()

print(len(data))

results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(results[:10])

print("END", time.time() - start)

start = time.time()

pool = mp.Pool(mp.cpu_count())

results = pool.starmap(howmany_within_range, [(row, 4, 8) for row in data])

pool.close()

print(len)

print("END2", time.time() - start)



