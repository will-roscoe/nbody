
import math
from multiprocessing import Pool, Process, Queue
from multiprocessing import freeze_support
def flatten(L):
    for item in L:
        try:
            yield from flatten(item)
        except TypeError:
            yield item
class Body:
    def __init__(self, a, b, c) -> None:
        self.a, self.b, self.c = a,b,c

def product(obj, q, i):
    q.put((f'{i}pr', obj.a+obj.c))
def sum(obj, q, i):
    q.put((f'{i}su', obj.a**obj.b))

bodies = flatten([[[Body(x,y,z) for x in range(5)] for y in range(5)] for z in range(5)])

if __name__ == "__main__":
    finished = []
    q = Queue()
    processes = []
    for i,b in enumerate(bodies):
        processes.append(Process(target=product,args=(b,q, i)))
        processes.append(Process(target=sum, args=(b,q, i)))
    for p in processes:
        p.start()
    for p in processes:
        finished.append(q.get(timeout=1))
    for p in processes:
        p.join()
    del processes, q

    for r in finished:
        print(r)