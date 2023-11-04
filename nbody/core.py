import math
import errors as e



def _vcomp(obj):
    if isinstance(obj, (HistoricVariable, HistoricVector, Vector)):
        return obj.c()
    else:
        return obj

NoneType = type(None)




#START of HistoricVariable Class
class HistoricVariable:
    def __init__(self, init_var, identity:str='HistoricVariable', units:str=None) -> None:
        if isinstance(init_var, (float, int)):
            self.hist = [init_var]
        elif isinstance(init_var, list):
            self.hist = init_var
        else: 
            e.raise_type_error('init_var', (int, float, list, tuple), init_var)
        if isinstance(identity, str):
            self.identity = identity
        else: 
            e.raise_value_error('identity', str, identity)
        if isinstance(units, str):
            self.units = units
        elif units is None:
            self.units = ''
        else: 
            e.raise_type_error('units', (NoneType, str), units)

    def c(self) -> float | int:
        return self.hist[-1]

    def next(self, next_val: int | float | list | tuple) -> None:
        if isinstance(next_val, (int, float)):
            self.hist.append(next_val)
        elif isinstance(next_val, (list, tuple)):
            for val in next_val:
                if isinstance(val, (float, int)):
                    self.hist.append(val)
                else: 
                    e.raise_list_type_error('next_val', (int, float), val)
        else: 
            e.raise_type_error('next_val', (int, float, list, tuple), next_val)

    # dUnder Methods
    def __len__(self) -> int:
        return len(self.hist)
    def __contains__(self, item:int | float) -> bool:
        return item in self.hist
    def __str__(self):
        return f'Historic Variable:(current={self.c()}, len={len(self)}) id="{self.identity}"'
    def __repr__(self) -> str:
        return f'HistoricVariable(current={self.c()} {self.units},\
len={len(self)}, hist={self.hist}, id={self.identity})'
    def __getitem__(self, ind: int | str) -> int | float | list | tuple:
        if isinstance(ind, int):
            return self.hist[ind]
        elif isinstance(ind, str):
            ind = ind.lower()
            if ind in ('first', 'initial'):
                return self.hist[0]
            elif ind in ('last', 'current'):
                return self.c()
            elif ind in ('full','hist', 'all'):
                return self.hist
            elif ind in ('past', 'old'):
                return self.hist[0:-1]
            else:
                e.raise_value_error('ind', str, ind)
        else:
            e.raise_type_error('ind', (str, int), ind)
    def __add__(self, other: int | float ) -> int | float:
        return self.c() + _vcomp(other)
    def __sub__(self, other: int | float) -> int | float:
        return self.c() - _vcomp(other)
    def __mul__(self, other: int | float) -> int | float:
        return self.c() * _vcomp(other)
    def __truediv__(self, other: int | float) -> int | float:
        return self.c() / _vcomp(other)
    def __floordiv__(self, other: int | float) -> int | float:
        return self.c() // _vcomp(other)
    def __iadd__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)):
            self.next(self.c() + temp)
            return self
        else:
            e.raise_type_error('other', (int, float), other)
    def __isub__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)):
            self.next(self.c() - temp)
            return self
        else:
            e.raise_type_error('other', (int, float), other)
    def __imul__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)):
            self.next(self.c() * temp)
            return self
        else:
            e.raise_type_error('other', (int, float), other)
    def __itruediv__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)):
            self.next(self.c() / temp)
            return self
        else:
            e.raise_type_error('other', (int, float), other)
#END OF HistoricVariable Class

#START of HistoricVector Class
class Vector:
    def __init__(self, li : tuple | list = None,
                x: float | int | list = None,
                y: float | int | list = None,
                z: float | int | list = None) -> None:
        if (isinstance(li, (tuple, list)) and len(li) == 3):
            self.X, self.Y, self.Z = li
        elif all(isinstance(var, (int, float)) for var in (x, y, z)):
            self.X, self.Y, self.Z = x, y, z
        else:
            e.raise_list_type_error('l,x,y,z', (tuple, list, int, float), (li,x,y,z))

    def c(self, usage:int =None) -> tuple:
        if usage is None:
            return (self.X, self.Y, self.Z)
        elif usage == 0:
            return self.X
        elif usage == 1:
            return self.Y
        elif usage == 2:
            return self.Z
        else:
            e.raise_out_of_range('c()', usage)

    def magnitude(self) -> float | int:
        return math.sqrt(sum([n**2 for n in self.c()]))
    def unit(self) -> tuple:
        if float(self.magnitude()) == 0.:
            return (0,0,0)
        else:
            return Vector(li=list(n/self.magnitude() for n in self.c()))
    def cross(self, other: tuple | list | object) -> tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            e1 = self.Y * temp[2] - self.Z * temp[1]
            e2 = self.Z * temp[0] - self.X * temp[2]
            e3 = self.X * temp[1] - self.Y * temp[0]
            return Vector(e1, e2, e3)
        else:
            e.raise_component_error('other', other)

    # dUnder Methods
    def __getitem__(self, ind: int | str) -> int | float | list | tuple:
        if isinstance(ind, str):
            ind = ind.lower()
            if ind in ('x', 'i'):
                return self.X
            elif ind in ('y', 'j'):
                return self.Y
            elif ind in ('z', 'k'):
                return self.Z
            elif ind == ('current'):
                return self.c()
            else:
                e.raise_value_error('ind', str, ind)
        else: 
            e.raise_type_error('ind', str, ind)
    def __str__(self) -> str: 
        return f'[{self.X}i+{self.Y}j+({self.Z})k]'
    def __repr__(self) -> str: 
        return f'{self.c()}'
    def __len__(self):
        return 1
    def __add__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            return Vector(li=[val + temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __sub__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            return Vector(li=[val - temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __mul__(self, other: list | tuple | float | int | object) -> list | tuple | float | int:
        temp = _vcomp(other)
        if isinstance(temp, (float, int)):
            return Vector(li=[val * temp for val in self['current']])
        elif len(temp) == 3:
            return sum(([val * temp[i] for i, val in enumerate(self['current'])]))
        else:
            e.raise_component_error('other or temp', temp)
    def __truediv__(self, other: int | float) -> list | tuple:
        temp = _vcomp(other)
        if isinstance(temp, (float, int)):
            return Vector(li=[val/temp for val in self['current']])
        else:
            e.raise_type_error('other', (float, int, object), other)

class HistoricVector:
    #init
    def __init__(self, x: float | int | list = None,
                 y: float | int | list = None,
                 z: float | int | list = None,
                 li : tuple | list = None, # update x from x_init etc
                 identity: str = None,
                 units_v: str = None) -> None:
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

        if isinstance(li, (tuple, list)) and len(li) == 3:
            self.X = HistoricVariable(li[0], f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(li[1], f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(li[2], f'{self.identity}_z', self.units)
        elif (isinstance(x, (int, float)) and
              isinstance(y, (int, float)) and
              isinstance(z, (int, float))):
            self.X = HistoricVariable(x, f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(y, f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(z, f'{self.identity}_z', self.units)
        else:
            e.raise_type_error('l,x,y,z', (tuple, list, float, int), (li,x,y,z))

    def x(self) -> float | int:
        return self.X.c()
    def y(self) -> float | int:
        return self.Y.c()
    def z(self) -> float | int:
        return self.Z.c()
    def c(self, usage:int =None) -> tuple:
        if usage is None:
            return (self.x(), self.y(), self.z())
        elif usage == 0:
            return self.x()
        elif usage == 1:
            return self.y()
        elif usage == 2:
            return self.z()
        else:
            e.raise_out_of_range('usage', usage)

    def magnitude(self) -> float | int:
        return math.sqrt(sum([n**2 for n in self.c()]))
    def unit(self) -> tuple:
        if float(self.magnitude()) == 0.:
            return (0,0,0)
        else: 
            return list(n/self.magnitude() for n in self.c())
    def cross(self, other: tuple | list | object) -> tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            e1 = self.y() * temp[2] - self.z() * temp[1]
            e2 = self.z() * temp[0] - self.x() * temp[2]
            e3 = self.x() * temp[1] - self.y() * temp[0]
            return (e1, e2, e3)
        else:
            e.raise_component_error('other',other)

    def next(self, next_vals: tuple | list) -> None:
        temp = _vcomp(next_vals)
        if isinstance(temp, (list, tuple)):
            if len(temp) == 3:
                self.X.next(temp[0])
                self.Y.next(temp[1])
                self.Z.next(temp[2])
            else:
                e.raise_component_error('next_vals', next_vals)
        else:
            e.raise_type_error('next_vals', (list, tuple), next_vals)

    # dUnder Methods
    def __len__(self):
        if len(self.X) == len(self.Y) and len(self.Y) == len(self.Z):
            return int(len(self.X))
        else:
            e.raise_unmatched_error(('x', 'y', 'z'),(self.X, self.Y, self.Z))
    def __str__(self) -> str:
        return f'"{self.identity}":[{self.x()}i+{self.y()}j+({self.z()})k], units ="{self.units}"'
    def __repr__(self) -> str:
        return f'HistoricVector("{self.identity}", current={self.c()}, lenx={len(self.X)})'
    def __getitem__(self, ind: int | str) -> int | float | list | tuple:
        if isinstance(ind, str):
            ind = ind.lower()
            if ind in ('x', 'i'):
                return self.x()
            elif ind in ('y', 'j'):
                return self.y()
            elif ind in ('z', 'k'):
                return self.z()
            elif ind == ('current'):
                return self.c()
            elif ind in ('full', 'all', 'hist'):
                return ((self.X.hist), (self.Y.hist), (self.Z.hist))
            elif ind in ('first', 'initial', 'past', 'old'):
                return (self.X[ind], self.Y[ind], self.Z[ind])
            else:
                e.raise_value_error('ind', str, ind)
        elif isinstance(ind, int):
            return (self.X[ind], self.Y[ind], self.Z[ind])
        else:
            e.raise_type_error('ind', (str, int), ind)
    def __add__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3: 
            return ([val + temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other',other)
    def __sub__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3: 
            return ([val - temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other',other)
    def __mul__(self, other: list | tuple | float | int | object) -> list | tuple | float | int:
        temp = _vcomp(other)
        if isinstance(temp, (float, int)): 
            return ([val * temp for val in self['current']])
        elif len(temp) == 3: 
            return sum(([val * temp[i] for i, val in enumerate(self['current'])]))
        elif isinstance(other, (list, tuple)): 
            e.raise_component_error('other',other)
        else:
            e.raise_type_error('other', (int, float, list, tuple), 'other')
    def __truediv__(self, other: int | float) -> list | tuple:
        temp = _vcomp(other)
        if isinstance(temp, (float, int)):
            return ([val/temp for val in self['current']])
        else:
            e.raise_type_error('other', (int, float), other)
#END of HistoricVector Class
