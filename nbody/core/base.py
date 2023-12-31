'''`nbody.base`
This library forms the internals of the structures and functioning of the rest of the objects users interact with.

Each `base.HistoricVector` instance is effectively comprised of 3 `base.HistoricVariable` instances, representing the
components of the vector. The classes have solely been built for use in a 3D instance,
but you may be able to use 2 dimensions inefficiently.

Arithmetic Functions for each of the object types all produce instances of `Variable` or `Vector` as the output,
except for magnitude.
'''

from __future__ import annotations
# ⤤ allows imports and references in a non chronological order.

#### Python Builtins
from numbers import Number
# ⤤ for type checking
#### 3rd Party Libs/Packages
from mpmath import mp, fp
# ⤤ for high prescision arithmetic
#### Local imports
from ..tools import errors as e
# ⤤ contains standardized error codes
from ..tools import _config as conf
# ⤤ MathContext object

math_conf = conf.MathContext(use=mp)

Iterable = (list, tuple)
Any = object
NoneType = type(None)
NumType = (Number, int, float ,type(mp.mpf(1)),type(fp.mpf(1)))
mp.pretty, fp.pretty = True, True


    
def _O(obj):  # noqa: N802
    if isinstance(obj,(Vector, Variable, HistoricVariable, HistoricVector)):
        return obj.c()
    else:
        return obj
    
def _V(obj):  # noqa: N802
    if isinstance(obj, Iterable):
        if len(obj) == 3:
            return Vector(li=obj)
        else:
            e.raise_len_error('obj in _V',3,obj)
    
    else:
        return obj

def _ntype(*objs):
    types = [type(_O(obj)) for obj in objs]
    for i,obj in enumerate(objs):
        if not isinstance(obj, NumType):
            return types[i]
    if all(ntype == int for ntype in types):
        return int
    else:
        return math_conf.type

    
def typecheck(argtype):
    if isinstance(argtype[0], Iterable):
        for arg in argtype:
            if not isinstance(arg[0],arg[1]):
                e.raise_type_error('input',arg[1],arg[0])
        
        return list(arg[0] for arg in argtype)
    else:
        if not isinstance(argtype[0],argtype[1]):
                e.raise_type_error('input',argtype[1],argtype[0])
        return argtype[0]




class Variable:
    def __init__(self,init_var,identity='Variable',units=None):
        (self.record,self.identity,self.units) = typecheck(((_ntype(init_var)(init_var),NumType),
                                                            (identity,str),(units,(str,NoneType))))
        self.units = (units if units else '')

    def c(self):
        return self.record
    
    def next(self, val):  # noqa: A003
        self.record = _O(val)

    def __len__(self):
        return 1
    
    def __contains__(self,item):
        return item in self.record

    def __str__(self):
        if self.identity in ('Variable', 'HistoricVariable'):
            return f'{self.c()} {self.units}, len:{len(self)}'
        else:
            return f'{self.c()} {self.units}, len:{len(self)} id:"{self.identity}"'
    
    def __repr__(self):
        if len(self) == 1:
            if self.identity in ('Variable', 'HistoricVariable'):
                return f'VarObj({self.c()} {self.units}, len={len(self)})'
            else:
                return f'VarObj({self.c()} {self.units}, len={len(self)}, id={self.identity})'
        elif self.identity in ('Variable', 'HistoricVariable'):
            return f'VarObj({self.c()} {self.units}, len={len(self)}, rec={self.record})'
        else:
            return f'VarObj({self.c()} {self.units}, len={len(self)}, rec={self.record}, id={self.identity})'
    
    def __add__(self,other):
        return Variable(math_conf.add(self.c(),_O(other)))
    
    __radd__ = __add__

    def __sub__(self,other):
        return Variable(math_conf.chop(math_conf.sub(self.c(),_O(other))))
    
    def __mul__(self,other):
        return Variable(math_conf.mul(self.c(),_O(other)))
    
    __rmul__ = __mul__

    def __truediv__(self,other):
        return Variable(math_conf.div(self.c(),_O(other)))
    
    def __iadd__(self,other):
        self.next(math_conf.add(self.c(),_O(other)))
        return self
    
    def __isub__(self,other):
        self.next(math_conf.chop(math_conf.sub(self.c(),_O(other))))
        return self
    
    def __imul__(self,other):
        self.next(math_conf.mul(self.c(),_O(other)))
        return self
    
    def __itruediv__(self,other):
        self.next(math_conf.div(self.c(),_O(other)))
        return self
    
    def __rsub__(self,other):
        return Variable(math_conf.chop(mp.fsub(_O(other), self.c())))
    
    
    def __rtruediv__(self,other):
        return Variable(math_conf.div(_O(other), self.c()))

    def __pow__(self, other):
        return Variable(math_conf.pow(self.c(),_O(other)))
    
    def __rpow__(self, other):
        return Variable(math_conf.pow(_O(other), self.c()))
    
    def __eq__(self, other):
        temp = _O(other)
        nt =_ntype(self.c(), temp)
        return str(nt(self.c())) == str(nt(temp)) or nt(self.c()) == nt(temp)
    
    def __lt__(self, other):
        temp = _O(other)
        nt =_ntype(self.c(), temp)
        return nt(self.c()) < nt(temp)
    
    
    def __le__(self, other):
        return (self.__eq__(other) or self.__lt__(other))
    
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    
    def __gt__(self, other):
        temp = _O(other)
        nt =_ntype(self.c(), temp)
        return nt(self.c()) > nt(temp)
    
    
    def __ge__(self, other):
        return (self.__eq__(other) or self.__gt__(other))
    
class HistoricVariable(Variable):
    def __init__(self,init_var,identity='HistoricVariable',units=''):
        if isinstance(init_var,NumType):
            self.record = [_ntype(init_var)(init_var)]
        
        elif isinstance(init_var,Iterable):
            self.record = [_ntype(*init_var)(val) for val in init_var]
        
        else: 
            e.raise_type_error('init_var',(*NumType,Iterable),init_var)
        
        self.type = type(self.record[0])
        (self.identity,self.units) = typecheck(((units,str),(identity,str)))

    
    def c(self):
        return self.record[-1]

    
    def next(self,next_val):  # noqa: A003
        temp = _O(next_val)
        if isinstance(temp,NumType):
            self.record.append(self.type(temp))
        
        elif isinstance(temp,Iterable):
            for val in temp:
                val = _O(val)
                if isinstance(val,NumType):
                    self.record.append(self.type(val))
                else: 
                    e.raise_list_type_error('next_val, temp',self.type,val)
        else: 
            e.raise_type_error('next_val, temp',(self.type,Iterable),next_val)

    
    def __len__(self):
        return len(self.record)
    
    
    def __getitem__(self,ind):
        
        if isinstance(ind,int):
            return self.record[ind]
        
        elif isinstance(ind,str):
            ind = ind.lower()
            if ind in ('first','initial'):
                return self.record[0]
            elif ind in ('last','current'):
                return self.c()
            elif ind in ('full','hist','all'):
                return self.record
            elif ind in ('past','old'):
                return self.record[0:-1]
            else:
                e.raise_value_error('ind',str,ind)
        
        else:
            e.raise_type_error('ind',(str,int),ind)








class Vector:
    
    def __init__(self,li=None,x=None,y=None,z=None):
        if li:
            (self.X,self.Y,self.Z) = [_ntype(*_O(li))(comp) for comp in _O(li)]
        elif all(isinstance(var, NumType) for var in (x,y,z)):
            [self.X,self.Y,self.Z] = [_ntype(x,y,z)(comp) for comp in (x,y,z)]
        else:
            e.raise_list_type_error('l,x,y,z',(Iterable,*NumType,*VectorType),(li,x,y,z))
    
    def c(self,usage=None):
        try:
            return {None:(self.X,self.Y,self.Z),0:self.X,1:self.Y,2:self.Z}[usage]
        except KeyError: 
            e.raise_out_of_range('c()',usage)

    def magnitude(self):
        return math_conf.norm(math_conf.matrix(self.c()))
    
    def unit(self):
        if float(self.magnitude()) == 0.:
            return NullVector()
        else:
            return Vector(self/self.magnitude())
    
    def __getitem__(self,ind):
        if isinstance(ind, str):
            try:
                return {'x':self.X,'i':self.X,'y':self.Y,'j':self.Y,'z':self.Z,'k':self.Z,'current':self.c()}[ind]
            except KeyError:
                e.raise_value_error('ind',str,ind)
        else: 
            e.raise_type_error('ind',str,ind)
    
    def __len__(self):
        return 3
    
    def __str__(self): 
        return f'{self.c()}'
    
    def __repr__(self): 
        return f'{self.c()}'
    
    def __iter__(self):
        return iter((self.X,self.Y,self.Z))
    def __add__(self,other):
        temp = _O(other)
        if len(temp) == 3:
            return Vector([math_conf.add(self.c(i),temp[i]) for i in range(3)])
        else:
            e.raise_component_error('other or temp',temp)
    
    __radd__ = __add__
    
    def __sub__(self,other):
        temp = _O(other)
        if len(temp) == 3:
            return Vector([math_conf.sub(self.c(i),temp[i]) for i in range(3)])
        else:
            e.raise_component_error('other or temp',temp)
    
    def __mul__(self,other):
        temp = _O(other)
        if isinstance(temp,Iterable) and len(temp) == 3:
            return Variable(math_conf.sum([math_conf.mul(val, temp[i]) for (i, val) in enumerate(self.c())]))
        elif isinstance(temp,NumType):
            _nt = math_conf.type
            return Vector([math_conf.mul(val,temp) for val in self.c()])
        else:
            e.raise_component_error('other or temp',temp)
    
    __rmul__ = __mul__

    def __truediv__(self,other):
        temp = _O(other)
        if isinstance(temp,NumType):
            return Vector([math_conf.div(val,temp) for val in self.c()])
        else:
            e.raise_type_error('other',NumType,other)

    
    def __rsub__(self,other):
        temp = _O(other)
        if len(temp) == 3:
            return Vector([math_conf.sub(temp[i], self.c(i)) for i in range(3)])
        else:
            e.raise_component_error('other or temp',temp)

    def __eq__(self, other):
        if other == None:
            return False
        else:
            temp = _O(other)
            nt = _ntype(*self.c(), *temp)
            [ex, ey, ez] =[str(nt(self.c(n))) == str(nt(temp[n])) or nt(self.c(n)) == nt(temp[n]) for n in range(3)]
            return ex and ey and ez
    


class HistoricVector(Vector):
    
    def __init__(self,x=None,y=None,z=None,li=None,identity=None,units_v=None):
        (self.identity,self.units) = typecheck(((identity,(NoneType,str)),(units_v,(NoneType,str))))
        self.units = (self.units if self.units is not None else '')
        self.identity =(self.identity if self.identity is not None else 'HistoricVector')
        li = ((x,y,z) if (isinstance(x,NumType) and isinstance(y,NumType) and isinstance(z,NumType)) else _O(li))
        nt = _ntype(*li)
        if isinstance(li,(tuple,list)) and len(li) == 3:
            (self.X,self.Y,self.Z) = list(HistoricVariable(vl,f'{self.identity}_{i}',self.units) for (vl,i) in 
                                          ((nt(li[0]),'x'),(nt(li[1]),'y'),(nt(li[2]),'z')))  
        else:
            e.raise_type_error('l,x,y,z',(Iterable,*NumType),(li,x,y,z))
        self.type = nt
    
    def x(self):
        return self.X.c()

    def y(self):
        return self.Y.c()

    def z(self):
        return self.Z.c()
    
    def c(self, usage=None):
        try:
            return {None:(self.x(),self.y(),self.z()),0:self.x(),1:self.y(),2:self.z()}[usage]
        except KeyError: 
            e.raise_out_of_range('c()',usage)
    
    def next(self,next_vals):  # noqa: A003
        temp = _O(next_vals)
        if isinstance(temp,Iterable):
    
            if len(temp) == 3:
                for i,comp in enumerate((self.X,self.Y,self.Z)):
                    comp.next(temp[i]) 
            else:
                e.raise_component_error('next_vals',next_vals)
        
        else:
            e.raise_type_error('next_vals',Iterable,next_vals)

    def __iter__(self):
        return iter((self.X.record, self.Y.record, self.Z.record))

    def __len__(self):
        if len(self.X) == len(self.Y) and len(self.Y) == len(self.Z):
            return int(len(self.X))
        else:
            e.raise_unmatched_error(('x', 'y', 'z'),(self.X, self.Y, self.Z))

    def __str__(self):
        return f'{self.c()}'
    
    def __repr__(self):
        return f'HistoricVector("{self.identity}", current={self.c()}, len={len(self)})'
    
    def __getitem__(self,ind):
        
        if isinstance(ind,str):
            ind = ind.lower()
            try:
                return {'current':self.c(),**dict.fromkeys(['x','i'],self.x()),**dict.fromkeys(['y','j'],self.y()), 
                           **dict.fromkeys(['z','k'],self.z()),**dict.fromkeys(['full','all','record'],
                                                                               ((self.X.record),(self.Y.record),(self.Z.record)))}[ind]
            except KeyError:
                if ind in ('first','initial','past','old'):
                    return (self.X[ind],self.Y[ind],self.Z[ind]) 
                else:
                    e.raise_value_error('ind',str,ind)
        
        elif isinstance(ind,int):
            return (self.X[ind],self.Y[ind],self.Z[ind])
        else:
            e.raise_type_error('ind',(str,int),ind)


class NullVector(Vector):
    def __init__(self):
        super().__init__(li=(0,0,0),x=None,y=None,z=None)

VectorType = (type(Vector((0,0,0))),type(HistoricVector(0,0,0)), type(NullVector()))
VarType = (type(Variable(0)),type(HistoricVariable(0)))