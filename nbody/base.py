from __future__ import annotations
# Python Builtins
import math
from decimal import Decimal

import numpy
# Local error definitions.
from . import errors as e

Any = object
NoneType = type(None)
DecType= type(Decimal())
NumType = (DecType, int, float)
Iterable = (list, tuple, numpy.ndarray)

def _O(obj):
    if isinstance(obj, (*VectorType, *VarType)):
        return obj.c()
    else:
        return obj
    
def _V(obj):
    if not isinstance(obj, VectorType):
        if len(obj) == 3:
            return Vector(li=obj)
        else:
            e.raise_len_error('obj in _V', 3, obj)
    else:
        return obj


def typecheck(argtype):
    for (inpt, typ) in argtype:
        if not isinstance(inpt, typ):
            e.raise_type_error('input', typ, inpt)

#START of Variable Class
class Variable:
    def __init__(self, init_var, identity='Variable', units=None):
        typecheck(((init_var, NumType), (identity, str), (units, (str, NoneType))))
        self.record = init_var
        self.identity = identity
        self.units = (units if units else '')
    def c(self) :
        return self.record
    def __len__(self):
        return 1
    def __contains__(self, item:Any) -> bool:
        return item in self.record
    def __str__(self):
        return f'{self.c()} {self.units}, len:{len(self)} id:"{self.identity}"'
    def __repr__(self):
        return f'VarType Object (current={self.c()} {self.units},\
len={len(self)}, rec={self.record}, id={self.identity})'
    def __add__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        return num_type(self.c()) + num_type(_O(other))
    def __sub__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        return num_type(self.c()) - num_type(_O(other))
    def __mul__(self, other):
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        return num_type(self.c()) * num_type(_O(other))
    def __truediv__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        return num_type(self.c()) / num_type(_O(other))
    def __floordiv__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        return num_type(self.c()) // num_type(_O(other))
    def __iadd__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        self.next(num_type(self.c()) + num_type(_O(other)))   
        return self
    def __isub__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        self.next(num_type(self.c()) - num_type(_O(other)))   
        return self
    def __imul__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        self.next(num_type(self.c()) * num_type(_O(other)))   
        return self
    def __itruediv__(self, other) :
        num_type = (type(_O(other)) if type(_O(other)) is not int else float)
        self.next(num_type(self.c()) / num_type(_O(other)))   
        return self
    
# --- SUBCLASSES ---
class HistoricVariable(Variable):
    def __init__(self, init_var, identity='HistoricVariable', units=None):
        if isinstance(init_var, NumType):
            self.record = [init_var]
        elif isinstance(init_var, list):
            self.record = init_var
        else: 
            e.raise_type_error('init_var', (*NumType, *Iterable), init_var)
        typecheck(((units, str), (identity, str)))
        self.identity = identity
        self.units = units

    def c(self) -> float | int:
        return self.record[-1]

    def next(self, next_val):
        if isinstance(next_val, NumType):
            self.record.append(next_val)
        elif isinstance(next_val, Iterable):
            for val in next_val:
                if isinstance(val, NumType):
                    self.record.append(val)
                else: 
                    e.raise_list_type_error('next_val', NumType, val)
        else: 
            e.raise_type_error('next_val', (*NumType, *Iterable), next_val)
    # dUnder Methods
    def __len__(self):
        return len(self.record)
    def __getitem__(self, ind):
        if isinstance(ind, int):
            return self.record[ind]
        elif isinstance(ind, str):
            ind = ind.lower()
            if ind in ('first', 'initial'):
                return self.record[0]
            elif ind in ('last', 'current'):
                return self.c()
            elif ind in ('full','hist', 'all'):
                return self.record
            elif ind in ('past', 'old'):
                return self.record[0:-1]
            else:
                e.raise_value_error('ind', str, ind)
        else:
            e.raise_type_error('ind', (str, int), ind)
#END OF Variable Class

#START of Vector Class
class Vector:
    def __init__(self, li = None,
                x = None,
                y = None,
                z = None):
        if li:
            li = _O(li)
        if (isinstance(li, Iterable) and len(li) == 3):
            self.X, self.Y, self.Z = li
        elif all(isinstance(var, NumType) for var in (x, y, z)):
            self.X, self.Y, self.Z = x, y, z
        else:
            e.raise_list_type_error('l,x,y,z', (*Iterable, *NumType, *VectorType), (li,x,y,z))

    def c(self, usage = None):
        _usage_lookup ={None: (self.X, self.Y, self.Z), 0:self.X, 1:self.Y, 2:self.Z}
        try:
            return _usage_lookup[usage]
        except KeyError: 
            e.raise_out_of_range('c()', usage)

    def magnitude(self):
        num_type = (type(self.c(0)) if type(self.c(0)) is not int else float)
        return num_type(math.sqrt(sum([n**2 for n in self.c()])))
    def unit(self):
        if float(self.magnitude()) == 0.:
            return Vector(li=(0,0,0))
        else:
            num_type =(type(self.c(0)) if type(self.c(0)) is not int else float)
            return Vector(li=list((num_type(n)/self.magnitude()) for n in self.c()))
    def cross(self, other):
        temp = _O(other)
        if len(temp) == 3:
            e1 = self.Y * temp[2] - self.Z * temp[1]
            e2 = self.Z * temp[0] - self.X * temp[2]
            e3 = self.X * temp[1] - self.Y * temp[0]
            return Vector(e1, e2, e3)
        else:
            e.raise_component_error('other', other)

    # dUnder Methods
    def __getitem__(self, ind):
        if isinstance(ind, str):
            _get_lookup = {'x':self.X, 'i':self.X, 'y':self.Y, 'j':self.Y, 'z':self.Z, 'k':self.Z, 'current':self.c()}
            try:
                return _get_lookup[ind]
            except KeyError:
                e.raise_value_error('ind', str, ind)
        else: 
            e.raise_type_error('ind', str, ind)
    def __len__(self):
        return 1
    def __str__(self): 
        return f'[{self.X}i+{self.Y}j+({self.Z})k]'
    def __repr__(self): 
        return f'{self.c()}, len={len(self)}'
    def __iter__(self):
        return iter((self.X, self.Y, self.Z))
    def __add__(self, other):
        temp = _O(other)
        if len(temp) == 3:
            num_type = (type(temp[0]) if type(temp[0]) is not int else float)
            return Vector(li=[num_type(val) + num_type(temp[i]) for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __sub__(self, other):
        temp = _O(other)
        if len(temp) == 3:
            num_type = (type(temp[0]) if type(temp[0]) is not int else float)
            return Vector(li=[num_type(val) - num_type(temp[i]) for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __mul__(self, other):
        temp = _O(other)
        if isinstance(temp, Iterable) and len(temp) == 3:
            num_type = (type(temp[0]) if type(temp[0]) is not int else float)
            return sum(([num_type(val) * num_type(temp[i]) for i, val in enumerate(self['current'])]))
        elif isinstance(temp, NumType):
            num_type = (type(temp) if type(temp) is not int else float)
            return Vector(li=[(num_type(val)*num_type(temp)) for val in self['current']])
        else:
            e.raise_component_error('other or temp', temp)
    def __truediv__(self, other):
        temp = _O(other)
        num_type = (type(temp) if type(temp) is not int else float)
        if isinstance(temp, NumType):
            return Vector(li=[(num_type(val)/num_type(temp)) for val in self['current']])
        else:
            e.raise_type_error('other', NumType, other)

# --- SUBCLASSES ---
class HistoricVector(Vector):
    def __init__(self, x = None,
                 y = None,
                 z = None,
                 li = None,
                 identity= None,
                 units_v= None):
        if isinstance(units_v, str):
            self.units = units_v
        elif units_v is None:
            self.units = ''
        else:
            e.raise_type_error('units_v', (str, NoneType), units_v)

        if isinstance(identity, str): 
            self.identity = identity
        elif identity is None: 
            self.identity = 'HistoricVector'
        else:
            e.raise_type_error('identity', (str, NoneType), identity)
        li = _O(li)
        if isinstance(li, (tuple, list)) and len(li) == 3:
            self.X = HistoricVariable(li[0], f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(li[1], f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(li[2], f'{self.identity}_z', self.units)
        elif (isinstance(x, NumType) and
              isinstance(y, NumType) and
              isinstance(z, NumType)):
            self.X = HistoricVariable(x, f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(y, f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(z, f'{self.identity}_z', self.units)
        else:
            e.raise_type_error('l,x,y,z', (*Iterable,*NumType), (li,x,y,z))

    def x(self) :
        return self.X.c()
    def y(self) :
        return self.Y.c()
    def z(self) :
        return self.Z.c()
        
    def c(self, usage=None):
        _usage_lookup ={None: (self.x(), self.y(), self.z()), 0:self.x(), 1:self.y(), 2:self.z()}
        try:
            return _usage_lookup[usage]
        except KeyError: 
            e.raise_out_of_range('c()', usage)
    def next(self, next_vals):
        temp = _O(next_vals)
        if isinstance(temp, Iterable):
            if len(temp) == 3:
                self.X.next(temp[0])
                self.Y.next(temp[1])
                self.Z.next(temp[2])
            else:
                e.raise_component_error('next_vals', next_vals)
        else:
            e.raise_type_error('next_vals', Iterable, next_vals)
    # dUnder Methods
    def __iter__(self):
        return iter((self.X.record, self.Y.record, self.Z.record))
    
    def __len__(self):
        if len(self.X) == len(self.Y) and len(self.Y) == len(self.Z):
            return int(len(self.X))
        else:
            e.raise_unmatched_error(('x', 'y', 'z'),(self.X, self.Y, self.Z))
    def __str__(self):
        
        return f'"{self.identity}":{self.x()}i+{self.y()}j+{self.z()})k, units ="{self.units}"'
    def __repr__(self):
        
        return f'HistoricVector("{self.identity}", current={self.c()}, len={len(self)})'
    def __getitem__(self, ind):
        if isinstance(ind, str):
            _get_lookup = {'current':self.c,**dict.fromkeys(['x', 'i'], self.x), **dict.fromkeys(['y', 'j'], self.y), 
                           **dict.fromkeys(['z', 'k'], self.z),**dict.fromkeys(['full', 'all', 'record'],
                                                                               ((self.X.record), (self.Y.record), (self.Z.record)))}
            ind = ind.lower()
            try:
                return _get_lookup[ind]()
            except KeyError:
                if ind in ('first', 'initial', 'past', 'old'):
                    return (self.X[ind], self.Y[ind], self.Z[ind])
                else:
                    e.raise_value_error('ind', str, ind)
        elif isinstance(ind, int):
            return (self.X[ind], self.Y[ind], self.Z[ind])
        else:
            e.raise_type_error('ind', (str, int), ind)
#END of Vector Class
VectorType = (Vector, HistoricVector)
VarType = (Variable, HistoricVariable)

#START of Body Class

            
            
            
            
            
            


#END of Simulation Class