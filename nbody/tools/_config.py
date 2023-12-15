from mpmath import fp, mp
from math import sqrt
def fltmat(obj):
    return obj
def fltnorm(*obj):
    return sqrt(sum(x**2 for x in obj))
class MathContext:
    def __init__(self,use=None):
        with open('MATHTYPE', 'r') as file:
            lines = file.readlines()
            use = eval(lines[0])
        print(f'{use}')
        if use == float:
            self.type = use       
            self.add = float.__add__
            self.sub = float.__sub__
            self.sum = sum
            self.div = float.__truediv__
            self.mul = float.__mul__
            self.pow = float.__pow__
            self.norm = fltnorm
            self.matrix = fltmat
            self.chop = float
        elif use in (mp, fp):
            self.type = use.mpf        
            self.add = use.fadd
            self.sum = use.fsum
            self.sub = use.fsub
            self.div = use.fdiv
            self.mul = use.fmul
            self.pow = use.power  
            self.norm = use.norm
            self.matrix = use.matrix
            self.chop = use.chop
        elif isinstance(use,dict):
            self.type = use['type']        
            self.add = use['add']
            self.sum = use['sum']
            self.sub = use['sub']
            self.div = use['div']
            self.mul = use['mul']
            self.pow = use['pow'] 
            self.norm = use['norm']
            self.matrix = use['matrix']
            self.chop = use['chop']
        else:
            raise TypeError
    def change(self, use=float):
        with open('MATHTYPE', 'x') as file:
            file.write(str(use))
        self.__init__(use)
                 