import multiprocessing as mproc
from random import choice, uniform, randint
objects = []
for i in range(10):
    ob = []
    for j in range(randint(1, 10)):
        ob.append(choice('abcdefghi'))
    objects.append(ob)



number = 0
#print(objects)


def compfunc(obj, number):
    
    slot0, slot1 = None, None
    if uniform(0.,1.)>0.7:
        return [obj[0:len(obj)//2], obj[len(obj)//2:]]
    else:
        return len(obj)


if __name__ == '__main__':
    pool = mproc.Pool(processes=4)  
    results = [pool.apply(compfunc, args=(obj, number)) for obj in objects]
    
    pool.close()
    pool.join()
    print(results)
    for result in results:
        if isinstance(result, list):
            [objects.append(res) for res in result]
        if isinstance(result, int):
            number += result

    print(number, objects)