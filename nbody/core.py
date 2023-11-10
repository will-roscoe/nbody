
import math

from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta
import astropy.units as u
import re
from decimal import Decimal

from tqdm import tqdm, trange
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from cycler import cycler

from . import errors as e

try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)

def _vcomp(obj):
    """
    Compute the component of an object, which can be a HistoricVariable, HistoricVector, or Vector.

    This function checks the type of the input object and returns the component if it's one of the specified types.
    
    Args:
        obj: An object that can be a HistoricVariable, HistoricVector, or Vector.

    Returns:
        The computed component of the input object if it matches one of the specified types. Otherwise, it returns the input object itself.
    """
    if isinstance(obj, (HistoricVariable, HistoricVector, Vector)):
        return obj.c()
    else:
        return obj


def sphere(pos, radius, N=20):
    (c, r) = (pos, radius)
    u, v = np.mgrid[0:2*np.pi:N*1j, 0:np.pi:N*1j]
    x = r*np.cos(u)*np.sin(v) + c[0]
    y = r*np.sin(u)*np.sin(v) + c[1]
    z = r*np.cos(v) + c[2]
    return x,y,z

NoneType = type(None)
DecType= type(Decimal('3.45'))
Numeric = (int, float, DecType)


#START of HistoricVariable Class
class HistoricVariable:
    """
    A class for managing historical numeric data with identity and units.

    This class allows you to store and manipulate historical numeric data, providing
    methods for accessing, updating, and performing mathematical operations on the data.

    Args:
        init_var (float, int, list, tuple): The initial value or list of values for the historic data.
        identity (str, optional): A descriptive label for the historic variable. Default is 'HistoricVariable'.
        units (str, optional): The units of the historic data. Default is None.

    Attributes:
        hist (list): A list containing the historical data.
        identity (str): A label for the historic variable.
        units (str): The units of the historic data.

    Methods:
        c() -> float | int: Get the current value of the historic data.
        next(next_val: int | float | list | tuple) -> None: Add the next value(s) to the historical data.
        __len__() -> int: Get the length of the historical data.
        __contains__(item: int | float) -> bool: Check if an item is present in the historical data.
        __str__() -> str: Return a string representation of the HistoricVariable.
        __repr__() -> str: Return a detailed string representation of the HistoricVariable.
        __getitem__(ind: int | str) -> int | float | list | tuple: Get specific elements from the historical data.
        __add__(other: int | float) -> int | float: Add the historic data to another value.
        __sub__(other: int | float) -> int | float: Subtract another value from the historic data.
        __mul__(other: int | float) -> int | float: Multiply the historic data by another value.
        __truediv__(other: int | float) -> int | float: Divide the historic data by another value.
        __floordiv__(other: int | float) -> int | float: Floor division of the historic data by another value.
        __iadd__(other: int | float): In-place addition of the historic data with another value.
        __isub__(other: int | float): In-place subtraction of the historic data with another value.
        __imul__(other: int | float): In-place multiplication of the historic data with another value.
        __itruediv__(other: int | float): In-place division of the historic data with another value.
    """
    def __init__(self, init_var, identity:str='HistoricVariable', units:str=None) -> None:
        if isinstance(init_var, Numeric):
            self.hist = [init_var]
        elif isinstance(init_var, list):
            self.hist = init_var
        else: 
            e.raise_type_error('init_var', (*Numeric, list, tuple), init_var)
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
        """
        Get the current value of the historic variable.

        Returns:
            float | int: The current value of the historic variable.
        """
        return self.hist[-1]

    def next(self, next_val: int | float | list | tuple) -> None:
        """
        Append the next value(s) to the historic variable.

        Args:
            next_val (int | float | list | tuple): The next value(s) to add to the historic variable.

        Raises:
            TypeError: If next_val is not of the expected type.
            ListTypeError: If an element in next_val is not of the expected type.
        """
        if isinstance(next_val, Numeric):
            self.hist.append(next_val)
        elif isinstance(next_val, (list, tuple)):
            for val in next_val:
                if isinstance(val, Numeric):
                    self.hist.append(val)
                else: 
                    e.raise_list_type_error('next_val', Numeric, val)
        else: 
            e.raise_type_error('next_val', (*Numeric, list, tuple), next_val)

    # dUnder Methods
    def __len__(self) -> int:
        """
        Get the length of the historic variable.

        Returns:
            int: The length of the historic variable.
        """
        return len(self.hist)
    def __contains__(self, item:int | float) -> bool:
        """
        Check if an item is present in the historic variable.

        Args:
            item (int | float): The item to check for.

        Returns:
            bool: True if the item is present, False otherwise.
        """
        return item in self.hist
    def __str__(self):
        """
        Get a string representation of the historic variable.

        Returns:
            str: A string representation of the historic variable.
        """
        return f'Historic Variable:(current={self.c()}, len={len(self)}) id="{self.identity}"'
    def __repr__(self) -> str:
        """
        Get a detailed string representation of the historic variable.

        Returns:
            str: A detailed string representation of the historic variable.
        """
        return f'HistoricVariable(current={self.c()} {self.units},\
len={len(self)}, hist={self.hist}, id={self.identity})'
    def __getitem__(self, ind: int | str) -> int | float | list | tuple:
        """
        Get the value(s) from the historic variable based on the index or label.

        Args:
            ind (int | str): The index or label to retrieve values.

        Returns:
            int | float | list | tuple: The value(s) from the historic variable.

        Raises:
            TypeError: If ind is not of the expected type.
            ValueError: If the provided string label (ind) is not recognized.
        """
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
        if isinstance(_vcomp(other), DecType) or isinstance(self.c(),DecType):
            return Decimal(self.c()) + Decimal(_vcomp(other))
        else:
            return self.c() + _vcomp(other)
    def __sub__(self, other: int | float) -> int | float:
        if isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            return Decimal(self.c()) - Decimal(_vcomp(other))
        else:
            return self.c() - _vcomp(other)
    def __mul__(self, other: int | float) -> int | float:
        if isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            return Decimal(self.c()) * Decimal(_vcomp(other))
        else:
            return self.c() * _vcomp(other)
    def __truediv__(self, other: int | float) -> int | float:
        if isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            return Decimal(self.c()) / Decimal(_vcomp(other))
        else:
            return self.c() / _vcomp(other)
    def __floordiv__(self, other: int | float) -> int | float:
        if isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            return Decimal(self.c()) // Decimal(_vcomp(other))
        else:
            return self.c() // _vcomp(other)
    def __iadd__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)) and isinstance(self.c(), (int, float)):
            self.next(self.c() + temp)
            return self
        elif isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            self.next(Decimal(self.c()) + Decimal(_vcomp(other)))   
            return self
        else:
            e.raise_type_error('other', Numeric, other)
    def __isub__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)) and isinstance(self.c(), (int, float)):
            self.next(self.c() - temp)
            return self
        elif isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            self.next(Decimal(self.c()) - Decimal(_vcomp(other)))   
            return self
        else:
            e.raise_type_error('other', Numeric, other)
    def __imul__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, (int, float)) and isinstance(self.c(), (int, float)):
            self.next(self.c() * temp)
            return self
        elif isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            self.next(Decimal(self.c()) * Decimal(_vcomp(other)))   
            return self
        else:
            e.raise_type_error('other', Numeric, other)
    def __itruediv__(self, other: int | float):
        temp = _vcomp(other)
        if isinstance(temp, Numeric):
            self.next(self.c() / temp)
            return self
        elif isinstance(_vcomp(other), DecType) or isinstance(self.c(), DecType):
            self.next(Decimal(self.c()) / Decimal(_vcomp(other)))   
            return self
        else:
            e.raise_type_error('other', Numeric, other)
#END OF HistoricVariable Class

#START of Vector Class
class Vector:
    """
    A class for representing three-dimensional vectors and performing vector operations.

    This class allows you to work with three-dimensional vectors, including operations
    for calculating vector magnitudes, unit vectors, and vector cross products.

    Args:
        li (tuple, list, optional): A tuple or list containing three components (X, Y, Z).
        x (float, int, list, optional): The X component of the vector.
        y (float, int, list, optional): The Y component of the vector.
        z (float, int, list, optional): The Z component of the vector.

    Attributes:
        X (float, int): The X component of the vector.
        Y (float, int): The Y component of the vector.
        Z (float, int): The Z component of the vector.

    Methods:
        c(usage: int = None) -> tuple | float | int: Get one or all components of the vector.
        magnitude() -> float | int: Calculate the magnitude of the vector.
        unit() -> tuple: Get the unit vector of the current vector.
        cross(other: tuple | list | object) -> tuple: Calculate the cross product of two vectors.

        Other dUnder Methods:
        __getitem__(ind: int | str) -> int | float | list | tuple: Get specific elements from the vector.
        __str__() -> str: Return a string representation of the Vector.
        __repr__() -> str: Return a detailed string representation of the Vector.
        __len__() -> int: Get the length of the Vector.
        __add__(other: list | tuple | object) -> list | tuple: Add another vector or components to the vector.
        __sub__(other: list | tuple | object) -> list | tuple: Subtract another vector or components from the vector.
        __mul__(other: list | tuple | float | int | object) -> list | tuple | float | int: Multiply the vector by another vector, scalar, or components.
        __truediv__(other: int | float) -> list | tuple: Divide the vector by a scalar.
    """
    def __init__(self, li : tuple | list = None,
                x: float | int | list = None,
                y: float | int | list = None,
                z: float | int | list = None) -> None:
        if (isinstance(li, (tuple, list)) and len(li) == 3):
            self.X, self.Y, self.Z = li
        elif all(isinstance(var, Numeric) for var in (x, y, z)):
            self.X, self.Y, self.Z = x, y, z
        elif isinstance(li, (Vector, HistoricVector)):
            (self.X, self.Y, self.Z) = li.c()
        else:
            print(type(li))
            e.raise_list_type_error('l,x,y,z', (tuple, list,*Numeric), (li,x,y,z))

    def c(self, usage:int =None) -> tuple:
        """
        Get the components of the vector.

        Args:
            usage (int, optional): The component to retrieve (0 for X, 1 for Y, 2 for Z). If None, returns all components.

        Returns:
            tuple: A tuple containing the X, Y, and Z components of the vector, or a single component if usage is specified.

        Raises:
            OutOfRangeError: If usage is not 0, 1, 2, or None.
        """
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
        """
        Calculate the magnitude (length) of the vector.

        Returns:
            float | int: The magnitude of the vector.
        """
        #try:
            #return math.sqrt(sum([n**2 for n in self.c()]))
        #except TypeError:
        if isinstance(self.c(0), DecType):
            return Decimal(math.sqrt(sum([n**2 for n in self.c()])))
        else:
            return math.sqrt(sum([n**2 for n in self.c()]))
    def unit(self) -> object:
        """
        Get the unit vector in the same direction as the current vector.

        Returns:
            Vector: A Vector representing the unit vector.
        """
        if float(self.magnitude()) == 0.:
            return Vector(li=(0,0,0))
        else:
            if isinstance(self.c(0), DecType):
                return Vector(li=list(Decimal(n)/Decimal(self.magnitude()) for n in self.c()))
            else:
                return Vector(li=list(n/self.magnitude() for n in self.c()))
    def cross(self, other: tuple | list | object) -> object:
        """
        Calculate the cross product of the vector with another vector or iterable.

        Args:
            other (tuple | list | object): Another vector or iterable to calculate the cross product with.

        Returns:
            Vector: The cross product of the vectors.

        Raises:
            ComponentError: If the other vector does not have three components.
        """
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
        """
        Get a component or components of the vector based on an index or label.

        Args:
            ind (int | str): The index or label to retrieve components (0 for X, 1 for Y, 2 for Z, or 'current').

        Returns:
            int | float | list | tuple: The component(s) of the vector.

        Raises:
            ValueError: If the provided string label (ind) is not recognized.
            TypeError: If ind is not of the expected type.
        """
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
        """
        Get a string representation of the vector.

        Returns:
            str: A string representation of the vector in the form '[Xi+Yj+(Z)k]'.
        """
        return f'[{self.X}i+{self.Y}j+({self.Z})k]'
    def __repr__(self) -> str: 
        """
        Get a detailed string representation of the vector.

        Returns:
            str: A detailed string representation of the vector as a tuple (X, Y, Z).
        """ 
        return f'{self.c()}'
    def __len__(self):
        return 1
    def __iter__(self):
        return iter((self.X, self.Y, self.Z))
    def __add__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            if isinstance(temp[0], DecType) or isinstance(self.c(0), DecType):
                return Vector(li=[Decimal(val) + Decimal(temp[i]) for i, val in enumerate(self['current'])])
            else:
                return Vector(li=[val + temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __sub__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3:
            if isinstance(temp[0], DecType) or isinstance(self.c(0), DecType):
                return Vector(li=[Decimal(val) - Decimal(temp[i]) for i, val in enumerate(self['current'])])
            else:
                return Vector(li=[val - temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other or temp', temp)
    def __mul__(self, other: list | tuple | float | int | object) -> list | tuple | float | int:
        temp = _vcomp(other)
        if isinstance(temp, (int, float)) and isinstance(self.c(0), (int, float)):
            return Vector(li=[val * temp for val in self['current']])
        if isinstance(temp, DecType) or isinstance(self.c(0), DecType):
            return Vector(li=[(Decimal(val)*Decimal(temp)) for val in self['current']])
        elif len(temp) == 3:
            if isinstance(temp[0], DecType) or isinstance(self.c(0), DecType):
                return sum(([Decimal(val) * Decimal(temp[i]) for i, val in enumerate(self['current'])]))
            else:
                return sum(([val * temp[i] for i, val in enumerate(self['current'])]))
        else:
            e.raise_component_error('other or temp', temp)
    def __truediv__(self, other: int | float) -> list | tuple:
        temp = _vcomp(other)
        if isinstance(temp, (int, float)) and isinstance(self.c(0), (int, float)):
            return Vector(li=[val/temp for val in self['current']])
        if isinstance(temp, DecType) or isinstance(self.c(0), DecType):
            return Vector(li=[(Decimal(val)/Decimal(temp)) for val in self['current']])
        else:
            e.raise_type_error('other', (*Numeric, object), other)
#END of Vector Class



#START of HistoricVector Class
class HistoricVector:
    """
    A class for managing historical three-dimensional vectors with identity and units.

    This class allows you to store and manipulate historical three-dimensional vectors
    composed of X, Y, and Z components. It provides methods for accessing individual
    components, performing vector operations, and managing historical data.

    Args:
        x (float, int, list, optional): The X component of the vector.
        y (float, int, list, optional): The Y component of the vector.
        z (float, int, list, optional): The Z component of the vector.
        li (tuple, list, optional): A tuple or list containing all three components (X, Y, Z).
        identity (str, optional): A descriptive label for the vector. Default is 'HistoricVector'.
        units_v (str, optional): The units of the vector components. Default is None.

    Attributes:
        X (HistoricVariable): A historic variable for the X component.
        Y (HistoricVariable): A historic variable for the Y component.
        Z (HistoricVariable): A historic variable for the Z component.
        identity (str): A label for the vector.
        units (str): The units of the vector components.

    Methods:
        x() -> float | int: Get the current X component of the vector.
        y() -> float | int: Get the current Y component of the vector.
        z() -> float | int: Get the current Z component of the vector.
        c(usage: int = None) -> tuple | float | int: Get one or all components of the vector.
        magnitude() -> float | int: Calculate the magnitude of the vector.
        unit() -> tuple: Get the unit vector of the current vector.
        cross(other: tuple | list | object) -> tuple: Calculate the cross product of two vectors.
        next(next_vals: tuple | list) -> None: Add the next values to the historical components.

        Other dUnder Methods:
        __len__() -> int: Get the length of the vector.
        __str__() -> str: Return a string representation of the HistoricVector.
        __repr__() -> str: Return a detailed string representation of the HistoricVector.
        __getitem__(ind: int | str) -> int | float | list | tuple: Get specific elements from the vector.
        __add__(other: list | tuple | object) -> list | tuple: Add another vector or components to the vector.
        __sub__(other: list | tuple | object) -> list | tuple: Subtract another vector or components from the vector.
        __mul__(other: list | tuple | float | int | object) -> list | tuple | float | int: Multiply the vector by another vector, scalar, or components.
        __truediv__(other: int | float) -> list | tuple: Divide the vector by a scalar.

    """
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
        elif (isinstance(x, Numeric) and
              isinstance(y, Numeric) and
              isinstance(z, Numeric)):
            self.X = HistoricVariable(x, f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(y, f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(z, f'{self.identity}_z', self.units)
        else:
            e.raise_type_error('l,x,y,z', (tuple, list,*Numeric), (li,x,y,z))

    def x(self) -> float | int:
        """
        Get the X component of the historic vector.

        Returns:
            float | int: The X component of the historic vector.
        """
        return self.X.c()
    def y(self) -> float | int:
        """
        Get the Y component of the historic vector.

        Returns:
            float | int: The Y component of the historic vector.
        """
        return self.Y.c()
    def z(self) -> float | int:
        """
        Get the Z component of the historic vector.

        Returns:
            float | int: The Z component of the historic vector.
        """
        return self.Z.c()
        
    def c(self, usage:int =None) -> tuple:
        """
        Get the components of the historic vector.

        Args:
            usage (int, optional): The component to retrieve (0 for X, 1 for Y, 2 for Z). If None, returns all components.

        Returns:
            tuple: A tuple containing the X, Y, and Z components of the historic vector, or a single component if usage is specified.

        Raises:
            OutOfRangeError: If usage is not 0, 1, 2, or None.
        """
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
        """
        Calculate the magnitude (length) of the historic vector.

        Returns:
            float | int: The magnitude of the historic vector.
        """
        return Decimal(math.sqrt(sum([n**2 for n in self.c()])))
    def unit(self) -> tuple:
        """
        Get the unit vector in the same direction as the historic vector.

        Returns:
            tuple: A tuple representing the unit vector.
        """
        if float(self.magnitude()) == 0.:
            return (0,0,0)
        else: 
            return list(Decimal(n)/self.magnitude() for n in self.c())
    def cross(self, other: tuple | list | object) -> tuple:
        """
        Calculate the cross product of the historic vector with another vector or iterable.

        Args:
            other (tuple | list | object): Another vector or iterable to calculate the cross product with.

        Returns:
            tuple: The cross product of the vectors.

        Raises:
            ComponentError: If the other vector does not have three components.
        """
        temp = _vcomp(other)
        if len(temp) == 3:
            e1 = self.y() * temp[2] - self.z() * temp[1]
            e2 = self.z() * temp[0] - self.x() * temp[2]
            e3 = self.x() * temp[1] - self.y() * temp[0]
            return (e1, e2, e3)
        else:
            e.raise_component_error('other',other)

    def next(self, next_vals: tuple | list) -> None:
        """
        Append the next values to the historic vector's components.

        Args:
            next_vals (tuple | list): The next values for the X, Y, and Z components of the historic vector.

        Raises:
            ComponentError: If next_vals does not have three components.
            TypeError: If next_vals is not of the expected type.
        """
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
    def __iter__(self):
        return iter((self.X.hist, self.Y.hist, self.Z.hist))
    def __len__(self):
        """
        Get the length of the historic vector.

        Returns:
            int: The length of the historic vector's components.

        Raises:
            UnmatchedError: If the X, Y, and Z components have different lengths.
        """
        if len(self.X) == len(self.Y) and len(self.Y) == len(self.Z):
            return int(len(self.X))
        else:
            e.raise_unmatched_error(('x', 'y', 'z'),(self.X, self.Y, self.Z))
    def __str__(self) -> str:
        """
        Get a string representation of the historic vector.

        Returns:
            str: A string representation of the historic vector in the form 'identity:[Xi+Yj+(Z)k], units ="units"'.
        """
        return f'"{self.identity}":[{self.x()}i+{self.y()}j+({self.z()})k], units ="{self.units}"'
    def __repr__(self) -> str:
        """
        Get a detailed string representation of the historic vector.

        Returns:
            str: A detailed string representation of the historic vector.
        """
        return f'HistoricVector("{self.identity}", current={self.c()}, lenx={len(self.X)})'
    def __getitem__(self, ind: int | str) -> int | float | list | tuple:
        """
        Get a component or components of the historic vector based on an index or label.

        Args:
            ind (int | str): The index or label to retrieve components (0 for X, 1 for Y, 2 for Z, or 'current').

        Returns:
            int | float | list | tuple: The component(s) of the historic vector.

        Raises:
            ValueError: If the provided string label (ind) is not recognized.
            TypeError: If ind is not of the expected type.
        """
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
            if isinstance(temp[0], DecType) or isinstance(self.x(), DecType):
                return ([Decimal(val) + Decimal(temp[i]) for i, val in enumerate(self['current'])])
            else:
                return ([val + temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other',other)
    def __sub__(self, other: list | tuple | object) -> list | tuple:
        temp = _vcomp(other)
        if len(temp) == 3: 
            if isinstance(temp[0], DecType) or isinstance(self.x(), DecType):
                return ([Decimal(val) - Decimal(temp[i]) for i, val in enumerate(self['current'])])
            else:
                return ([val - temp[i] for i, val in enumerate(self['current'])])
        else:
            e.raise_component_error('other',other)
    def __mul__(self, other: list | tuple | float | int | object) -> list | tuple | float | int:
        temp = _vcomp(other)
        if isinstance(temp, (int,float)) and isinstance(self.x(), (int, float)): 
            return ([val * temp for val in self['current']])
        if isinstance(temp, DecType) or isinstance(self.x(), DecType):
            return ([(Decimal(val) * Decimal(temp)) for val in self['current']])
        elif len(temp) == 3: 
            if isinstance(temp[0], DecType) or isinstance(self.x(), DecType):
                return sum(([Decimal(val) * Decimal(temp[i]) for i, val in enumerate(self['current'])]))
            else:
                return sum(([val * temp[i] for i, val in enumerate(self['current'])]))
        elif isinstance(other, (list, tuple)): 
            e.raise_component_error('other',other)
        else:
            e.raise_type_error('other', (int, float, list, tuple), 'other')
    def __truediv__(self, other: int | float) -> list | tuple:
        temp = _vcomp(other)
        if isinstance(temp, (int,float)) and isinstance(self.x(), (int, float)): 
            return ([val/temp for val in self['current']])
        elif isinstance(temp, DecType) or isinstance(self.x(), DecType):
            return ([(Decimal(val)/Decimal(temp)) for val in self['current']])
        else:
            e.raise_type_error('other', Numeric, other)
#END of HistoricVector Class

#START of Body Class
class Body:
    """
    A class representing a physical body with mass, position, velocity, and acceleration.

    This class allows you to define a physical body with attributes such as mass, position,
    velocity, and acceleration. It provides methods to update and evaluate the body's
    state over time.

    Args:
        mass (float, int): The mass of the body (in kilograms).
        init_pos (list, tuple): The initial position of the body (in meters).
        init_vel (list, tuple, optional): The initial velocity of the body (in meters per second).
        radius (float, int, optional): The radius of the body (in meters).
        identity (str, optional): A descriptive label for the body. Default is auto-generated.

    Attributes:
        mass (HistoricVariable): A historic variable for the mass of the body (in kilograms).
        pos (HistoricVector): A historic vector for the position of the body (in meters).
        vel (HistoricVector): A historic vector for the velocity of the body (in meters per second).
        acc (HistoricVector): A historic vector for the acceleration of the body (in meters per second squared).
        identity (str): A label for the body.

    Methods:
        evaluate(dt: int | float = 1) -> None: Evaluate the body's state over a time interval (dt).
        update(dt: int | float = 1, vel_change: list | tuple = None, acc_change: list | tuple = None) -> None: Update the body's state based on changes in velocity and acceleration.

        Other dUnder Methods:
        __str__() -> str: Return a string representation of the Body.
        __repr__() -> str: Return a detailed string representation of the Body.

    """
    def __init__(self,
                mass: float | int,
                init_pos: list | tuple,
                init_vel: list | tuple=(0,0,0),
                radius: float | int = 0,
                bounce: float | int = 0.999,
                color: str= None,
                identity:str = None) -> None:
        if isinstance(identity, str):
            self.identity = identity
        elif identity is None:
            self.identity = f'body[m={mass}]'
        else:
            e.raise_type_error('identity', str, identity)
        self.acc = HistoricVector(0,0,0,
                             identity=f'{self.identity}_acc',
                             units_v='ms^-2')
        if isinstance(init_pos, (list, tuple)):
            self.pos = HistoricVector(li=init_pos,
                                 identity=f'{self.identity}_pos',
                                 units_v='m')
        else:
            e.raise_type_error('init_pos', (list, tuple), init_pos)

        if isinstance(init_vel, (list, tuple)):
            self.vel = HistoricVector(li=init_vel,
                                 identity=f'{self.identity}_vel',
                                 units_v='ms^-1')
        else:
            e.raise_type_error('init_vel', (list, tuple), init_vel)
        if isinstance(mass, Numeric):
            self.mass = HistoricVariable(mass,
                                identity=f'{self.identity}_mass',
                                units='kg')
        else:
            e.raise_type_error('mass', Numeric, mass)
        if isinstance(radius, Numeric):
            self.radius = HistoricVariable(radius,
                                identity=f'{self.identity}_rad',
                                units='m')
        else:
            e.raise_type_error('radius', Numeric, radius)
        if isinstance(bounce, Numeric):
            self.bounce = bounce
        else:
            e.raise_type_error('bounce', Numeric, bounce)
        if isinstance(color, (NoneType, str)):
            self.color = color
        else:
            e.raise_type_error('color',(NoneType, str), color)

    def __str__(self) -> str:
        return f'Body("{self.identity}",\n\
    mass="{self.mass.c()} {self.mass.units}",\n\
    currentpos="{self.pos.c()} {self.pos.units}",\n\
    currentvel="{self.vel.c()} {self.vel.units}",\n\
    currentvel="{self.acc.c()} {self.acc.units}")'
    def __repr__(self) -> str:
        return f'Body("{self.identity}", m={self.mass.c()}, r={self.pos.c()},\
v={self.vel.c()}), a={self.acc.c()})'


    def evaluate(self, dt: int | float =1) -> None:
        """
        This method updates the position, velocity, and acceleration of a body over
        a small time interval `dt`. It ensures that the lengths of position and
        velocity vectors are consistent and calculates the acceleration.

        Args:
            dt (int or float, optional): The time step for the evaluation (default is 1).

        Returns:
            None

        Note:
            This method modifies the internal state of the body object.
        """
        while len(self.pos) - len(self.vel) != 0:
            len_check = len(self.pos) - len(self.vel)
            if len_check < 0:
                for i in range(len_check, 0, 1):
                    self.pos.next(self.pos + Vector((self.vel[i]))*dt)
            elif len_check > 0:
                for i in range (-len_check, 0, 1):
                    self.vel.next(Vector((self.pos[len_check]-self.pos[len_check-1]))/dt)
        while len(self.pos) > len(self.acc):
            self.acc.next((0,0,0))


    def update(self, dt: int | float =1,
                vel_change: list | tuple=None,
                acc_change: list | tuple=None,
                vel_next: list | tuple=None,
                acc_next: list | tuple=None, ) -> None:
        """
        Update the position and velocity of a body over a small time interval.

        This method updates the position and velocity of a body based on the time step `dt`, 
        velocity changes (`vel_change`), and acceleration changes (`acc_change`). It ensures 
        that the lengths of position and velocity vectors are consistent.

        Args:
            dt (int or float, optional): The time step for the update (default is 1).
            vel_change (list or tuple, optional): Velocity changes as a 3-element list or tuple 
            (default is None).acc_change (list or tuple, optional): Acceleration changes as
            a 3-element list or tuple (default is None).

        Returns:
            None

        Raises:
            ComponentError: If `vel_change` or `acc_change` does not have 3 elements.
            TypeError: If `vel_change` or `acc_change` is not a list or tuple.
            EvaluationError: If the length of the position and velocity vectors is inconsistent.

        Note:
            This method modifies the internal state of the body object.

        Examples:
            # Create a body object
            body = Body()

            # Update the body's position and velocity with a time step of 0.1 seconds and velocity change
            body.update(0.1, vel_change=(1.0, 0.0, 0.0))

            # Update the body's position and velocity with a time step of 0.1 seconds and acceleration change
            body.update(0.1, acc_change=(0.0, -9.81, 0.0))
        """
        if vel_next:
            vel = vel_next
        else:
            vel = self.vel.c()
        if acc_next:
            acc = acc_next
        else:
            acc = self.acc.c()
            
        if acc_change is not None:
            if isinstance(acc_change, (list, tuple)):
                if len(acc_change) == 3:
                    self.vel.next((Vector((acc_change))*dt + vel).c())
                    self.acc.next(acc_change)
        else:
            self.vel.next((Vector((self.acc*dt))+vel).c())
            self.acc.next((0,0,0))
        if vel_change is not None:
            if isinstance(vel_change, (list, tuple)):
                if len(vel_change) == 3:
                    self.pos.next(self.pos + (Vector(vel) + vel_change)*dt)
        else:
            self.pos.next(self.pos + Vector(vel)*dt)
        #self.evaluate(dt)
    
    def _reinitialise(self, init_pos=None, init_vel=None):
        self.acc = HistoricVector(0,0,0,
                             identity=f'{self.identity}_acc',
                             units_v='ms^-2')
        if init_pos != None:
            if isinstance(init_pos, (list, tuple)):
                self.pos = HistoricVector(li=init_pos,
                                    identity=f'{self.identity}_pos',
                                    units_v='m')
            else:
                e.raise_type_error('init_pos', (list, tuple), init_pos)
        if init_vel != None:
            if isinstance(init_vel, (list, tuple)):
                self.vel = HistoricVector(li=init_vel,
                                    identity=f'{self.identity}_vel',
                                    units_v='ms^-1')
            else:
                e.raise_type_error('init_vel', (list, tuple), init_vel)
#END of Body Class


def horizons_query(searchquery: str,
                   observer: str = '0',
                   time: str = '2023-11-03',
                   num_type: type = float,
                   return_type: str = 'body'):
    """
    Create a Horizons object representing a celestial body using Horizons query data.

    This function queries NASA's Horizons system for information about a celestial body and creates a
    `Body` object with the retrieved data, including mass, initial position, initial velocity, radius,
    and identity.

    Args:
        searchquery (str): The identifier or name of the celestial body to query.
        observer (str, optional): The observer location (default is '0' for the solar system barycenter).
        time (str, optional): The date for which to retrieve data in the 'YYYY-MM-DD' format
                                (default is '2023-11-03').
        num_type (type, optional): Numeric Type to format values from source as. (default is float)
        return_type (str, optional): what format to return the object in. either 'print' to print output, 
                                    'dict' to return dictionary containing values or 'body' to return a Body
                                    instance. (default is 'body')
    Returns:
        Body: A `Body` object representing the celestial body with mass, initial position, initial velocity,
        radius, and identity information.
    OR  dict: dictionary with keys ('identity', 'mass', 'radius', 'position', 'velocity'). last 2 are 3-tuples.
    OR  print: prints a string representation of the dict above in a user readable form.
    Raises:
        TypeError: If one or more args are incorrect type.
        ValueError: if 'return_type' is not a valid option. make sure it is any of 'print', 'dict', 'body'.

    """
    _args = ((searchquery, str),(observer,str), (time, str), (num_type, type), (return_type,str))
    for (inpt, typ) in _args:
        if not isinstance(inpt, typ):
            e.raise_type_error('input', typ, inpt)
    if all([return_type.lower() != opt for opt in ('body', 'dict', 'print')]) :
        e.raise_value_error('return_type', type, return_type)
    
    tqdm.write(f'Querying "{searchquery}" @JPL Horizons System')
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
    _raw = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors(get_raw_response=True)
    _tab = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors()
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _raw)
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        tqdm.write(f'could not find name for object "{searchquery}", reverting to placeholder name.')
        name = f'JPL_{searchquery}'
    
    _raw =_raw.split('Ephemeris')[0]

    m_exp_unit = _raw.split('Mass')[1].split('^')[1].split('=')[0].strip(') ,~ (')
    m_exp = "".join(c for c in m_exp_unit if c.isdigit() or c == '.')
    m_unit = "".join(c for c in m_exp_unit.lower() if c == 'k' or c == 'g')
    mass = "".join(c for c in  _raw.split('Mass')[1].split('^')[1].split('=')[1].strip(') ,~ (').lower() if 
                   c.isdigit() or c == '.' or c=='+' or c=='-')
    
    try:
        rad_string = _raw.lower().split('vol. mean radius')[1].split('=')[0:2]
    except IndexError:
        rad_string = _raw.lower().split('radius')[1].split('=')[0:2] 
    
    r_unit = "".join(c for c in rad_string[0].lower() if c == 'k' or c == 'm') 
    
    rad = list(c for c in rad_string[1].split(' ') if c != '') 
    _rad = []
    for x in rad:
        if any(char.isdigit() for char in x): 
            _rad.append(x)
    
    mass = num_type(''.join(mass.split('+-')[0]))*num_type(10**int(m_exp))
    rad = num_type(_rad[0].split('+-')[0])

    if r_unit == 'km':
        rad *= 1000
    if m_unit == 'g':
        mass /= 1000

    x, y, z = [num_type(_tab[pos].quantity[0].to_value(u.m)) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [num_type(_tab[pos].quantity[0].to_value(u.m/u.s)) for pos in ('vx', 'vy', 'vz')]
    

    
    if return_type.lower() == 'print':
        print(f'''Object: {name}\n***********\nMass: {mass} kg\nRadius: {rad} m\n
***********\nposition: ({x}, {y}, {z}) (m)\nvelocity: ({vx}, {vy}, {vz}) (ms^-1)\n
***********\nQuery Date: {time}''')
    elif return_type.lower() == 'dict':
        return {'identity': name, 'mass': mass, 'radius': rad, 'position':(x,y,z), 'velocity':(vx,vy,vz)}
    elif return_type.lower() == 'body':
        return Body(mass=mass, init_pos=(x,y,z), init_vel=(vx,vy,vz), radius=rad, identity=name)


def horizons_batch(search_queries: tuple | list,
                   observer: str = '0',
                   time: str = '2023-11-03',
                   num_type: type = float,
                   return_type: str = 'body'):
    new_bodies = []
    for query in tqdm(search_queries, desc='Getting data from JPL Horizons', unit='queries'):
        if isinstance(query, str):
            new_bodies.append(horizons_query(query, observer, time, num_type, return_type))
        else:
            raise e.raise_type_error('item in search_queries', str, query)
    return new_bodies

#START of PhysEngine class
class PhysEngine:
    """
    A physics engine for simulating the dynamics of bodies in a physical environment.

    This class represents a physics engine that can simulate the behavior of multiple
    bodies with their interactions, gravity, and collisions. It allows you to attach
    bodies and create planes to simulate a physical environment.

    Args:
        dt (int, float): The time step for the simulation (in seconds).

    Attributes:
        bodies (list): A list of `Body` objects representing the physical bodies in the simulation.
        dt (int, float): The time step for the simulation (in seconds).
        planes (list): A list of planes that can be used for collision detection.

    Methods:
        attach_bodies(new_bodies): Attach new bodies to the simulation.
        create_plane(const_axis='z', const_val=0): Create a plane for collision detection.
        check_collision(body): Check for collisions between a body and other bodies or planes.
        find_gravity(): Calculate the gravitational forces between bodies.
        evaluate(): Perform the physics simulation for the current time step.

    """
    def __init__(self, dt: int | float = 1, checking_range:int=3):
        self.bodies = []
        if isinstance(dt, (Numeric)):
            self.dt = dt
        self.planes = []
        self.fields = []
        self._rangechk = checking_range
    def attach_bodies(self, new_bodies):
        """
        Attach a list of Body objects to the physics engine.

        Args:
            new_bodies (list or tuple): A list of Body objects to attach.

        Raises:
            TypeError: If new_bodies is not a list or tuple.
            TypeError: If an element in new_bodies is not of type Body.
        """
        if isinstance(new_bodies, (list, tuple)):
            for i, new_body in enumerate(new_bodies):
                if isinstance(new_body, Body) or issubclass(type(new_body), Body):
                    self.bodies.append(new_body)
                else:
                    e.raise_type_error(f'new_body at index {i}', type(Body), new_body)
        else:
            e.raise_type_error('new_bodies', (list, tuple), new_bodies)
        tqdm.write(f'{len(self.bodies)} bodies attached.')
    def make_relative_to(self, target_body):
        for body in self.bodies:
            if body != target_body:
                body._reinitialise((body.pos-target_body.pos), (body.vel-target_body.vel))
                body.identity = f'{body.identity}:Rel.{target_body.identity}'
        target_body._reinitialise((0,0,0), (0,0,0))
        tqdm.write(f"Bodies' positions and velocities have been made relative to {target_body.identity}.")
        target_body.identity = f'{target_body.identity}(Static)'
  
    def orbit_around(self, main_body, other_bodies):
            pass
    def create_acceleration(self, accel_vector):
        if len(accel_vector) == 3 and isinstance(accel_vector[0], Numeric):
            self.fields+= [Vector(li=accel_vector)] 
    def create_plane(self, const_axis='z', const_val = 0):
        """
        Create a collision plane in the physics engine.

        Args:
            const_axis (str, optional): The constant axis for the collision plane ('x', 'y', or 'z').
            const_val (int or float, optional): The constant value for the collision plane.

        Raises:
            ValueError: If const_axis is not 'x', 'y', or 'z', or if const_val is not an int or float.
        """
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, Numeric):
            self.planes.append([const_axis, const_val])
            tqdm.write(f'constant plane {const_axis}={const_val} has been initialized.')
        else:
            e.raise_value_error('const_axis,const_val',((str),(int,float)),(const_axis,const_val))
    def _check_collision(self, body, co_restitution=0):
        """
        Check for collisions between a body and other bodies or planes.

        Args:
            body (Body): The body to check for collisions.

        Returns:
            tuple: The resulting change in velocity for the body after collisions.

        Note:
            This method calculates and handles collisions with other bodies and collision planes.
        """
        returned_coll = False
        range_chk = 3
        for bod in self.bodies:
            if body != bod: 
                body_dists = list(((Vector(bod.pos.c()) + Vector(bod.vel.c())*m*self.dt) - 
                                (Vector(body.pos.c()) + Vector(body.vel.c())*m*self.dt)).magnitude() 
                                for m in range(1))
                if any([(body_dists[i] < (bod.radius.c() + body.radius.c())) for i in range(1)]):
                    returned_coll = True
                    n = Vector(bod.pos-body.pos)/body_dists[0]
                    meff = 1/((1/bod.mass.c())+(1/body.mass.c()))
                    vimp = n*(body.vel - bod.vel)
                    imp = vimp*(1+co_restitution)*meff
                    dv = n*(-imp/body.mass.c())
                    return (dv+body.vel), False
        
        unitcomps = {'x':Vector((1,0,0)), 'y':Vector((0,1,0)), 'z':Vector((0,0,1))}
        for [pl_ax, pl_val] in self.planes:
            body_est_dists =  list(abs((Vector(body.pos.c()) +
            Vector(body.vel.c())*m*self.dt)[pl_ax] - pl_val) for m in range(self._rangechk))
            body_cur_dist = abs(body.pos[pl_ax] - pl_val)
            if any([(body_est_dists[i] < body.radius.c()) for i in range(self._rangechk)]) or\
                body_cur_dist < body.radius.c():
                if  body_cur_dist/body.radius.c() <= 1.01 and body_est_dists[-1]/body.radius.c() <= 1.01:
                    on_plane = True
                else:
                    on_plane = False
                pl_norm = unitcomps[pl_ax]
                returned_coll = True
                return (pl_norm*(-2*(body.vel*pl_norm)) + body.vel)*co_restitution, on_plane

        if not returned_coll:
            return Vector(body.vel.c()), False
    

    def _find_gravity(self, body):
        """
        Calculate and apply gravitational forces between bodies.

        Note:
            This method calculates gravitational forces and updates the acceleration of each body.
        """
        bod1_data = [body, HistoricVector(0,0,0)]
        for bod2 in self.bodies: 
            if bod2 != body:
                temp_dist = body.pos - bod2.pos
                dist_12 = Vector(temp_dist)
                unit_12 = dist_12.unit()
                if isinstance(body.mass.c(), DecType) or isinstance(bod2.mass.c(), DecType):
                    b1m, b2m = Decimal(body.mass.c()), Decimal(bod2.mass.c())
                    force_on_bod1_t = (-1*Decimal(G) * b1m * b2m)
                    force_on_bod1 = unit_12*(force_on_bod1_t/Decimal(dist_12.magnitude())**2)
                else:
                    force_on_bod1_t = (-1*G * body.mass.c() * bod2.mass.c())
                    force_on_bod1 = unit_12*(force_on_bod1_t/(dist_12.magnitude())**2)
                bod1_data[1] = HistoricVector(li=(bod1_data[1] + force_on_bod1))
        return bod1_data[1]/body.mass.c()

    def evaluate(self):
        """
        Evaluate the motion of attached bodies in the physics engine.

        Note:
            This method calculates the motion of attached bodies based on gravitational forces and collisions.
        """
        _temp = [[0,0,0] for _ in self.bodies]
        for i, body in enumerate(self.bodies):
            _temp[i][0:2] = self._check_collision(body, body.bounce)
        for i, body in enumerate(self.bodies):
            _temp[i][2] = self._find_gravity(body)
        for i, body in enumerate(self.bodies):
            col_vel, on_plane, acc_g = _temp[i]
            if not on_plane:
                fieldvel = list(sum(f.c(i) for f in self.fields) for i in range(3))
            else:
                fieldvel = (0,0,0)
            body.update(self.dt, vel_next=(col_vel+fieldvel).c(), acc_change=acc_g)
#END of PhysEngine class
from multiprocessing import Pool
class PhysEngineMP(PhysEngine):
    def __init__(self, dt: int | float = 1, checking_range: int = 3):
        super().__init__(dt, checking_range)

    def _process_body(self,body: Body):
        return [body, *self._check_collision(body), self._find_gravity(body)]


    def evaluate(self):
        if __name__ == '__main__':
            with Pool as workers:
                _temp = workers.map(self._process_body, self.bodies)
        for _entry in _temp:
            body ,col_vel, on_plane, acc_g = _entry
            if not on_plane:
                fieldvel = list(sum(f.c(i) for f in self.fields) for i in range(3))
            else:
                fieldvel = (0,0,0)
            body.update(self.dt, vel_next=(col_vel+fieldvel).c(), acc_change=acc_g)
# START of Simulation Class
class Simulation:
    """
    A class for visualizing and controlling a physics simulation.

    This class allows you to create a 3D visualization of a physics simulation using
    Matplotlib and provides controls for adjusting the simulation parameters and viewing
    options. It includes features such as zoom, grid visibility, shadows, acceleration and
    velocity vectors, body models, and labeling options.

    Args:
        name (str): The name of the simulation.
        engine (PhysEngine|NoneType): The physics engine to simulate the physics of the system.
        focus_body (Body|NoneType): The body to focus on during the simulation.
        focus_range (int, float): The range for focusing on the body.
        autoscale (bool): Enable or disable autoscaling.
        show_grid (bool): Show or hide the grid in the visualization.
        show_shadows (bool): Show or hide shadows.
        show_acceleration (bool): Show or hide acceleration vectors.
        show_velocity (bool): Show or hide velocity vectors.
        vector_size (int, float): The size of velocity and acceleration vectors.
        labelling_type (str): The type of labeling, either 'legend' or 'label'.
        body_model (str): The model for representing bodies, such as 'dots', 'wireframe', or 'surface'.
        guistyle (str): The GUI style, either 'default' or 'dark'.

    Methods:
        start(frames, interval, duration, fps): Start the simulation animation.
    
    """
    def __init__(self,
                name: str = 'Nbody Simulation',
                engine: PhysEngine|NoneType = None,
                focus_body: Body|NoneType = None,
                focus_range: int|float|NoneType = None,
                autoscale: bool = True,
                show_grid: bool = True, 
                show_shadows: bool = False,
                show_acceleration: bool = False,
                show_velocity: bool = False,
                vector_size: int|float = 1,
                labelling_type: str = 'legend',
                body_model: str = 'dots',
                guistyle: str = 'default'):
        _argspec = {engine:(PhysEngine),focus_range:(*Numeric, NoneType),
autoscale:bool,show_grid:bool,show_shadows:bool,show_acceleration:bool,show_velocity:bool,
vector_size:Numeric,labelling_type:str,body_model:str,guistyle:str}
        
            
        for i,(arg,typ) in enumerate(_argspec.items()):
            if not isinstance(arg, typ):
                e.raise_type_error(f'Simulation arg({i})', typ, arg)
                break
        else:
            self._engine = engine
            self.focus_range = focus_range
            self.autoscale = autoscale
            self.show_grid = show_grid
            self.show_shadows = show_shadows
            self.show_acceleration = show_acceleration
            self.show_velocity = show_velocity
            self.vector_size = vector_size
            self.labelling_type = labelling_type
            self.body_model = body_model
            self.guistyle = guistyle
        if issubclass(type(focus_body), Body) or isinstance(focus_body, (NoneType, Body)):
            self.focus_body = focus_body
        else:
            e.raise_type_error('focus_body', (Body, 'subclassBody', NoneType), focus_body)

        
        self.fig = plt.figure(name, figsize=(16,9))
        self.ax = self.fig.add_subplot(projection = '3d')
        self.ax.computed_zorder = False
        self.ax1 = self.fig.add_axes((0.05,0.25,0.05,0.5))
        self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.zoom_slider = Slider(self.ax1, 'Zoom', 0.1, 10, valinit=1, orientation='vertical')
        self._frameskip, self._plotskip = 1, 1

        def _clearpanes():
            self.ax.xaxis.set_pane_color((0.,0.,0.,0.))
            self.ax.yaxis.set_pane_color((0.,0.,0.,0.))
            self.ax.zaxis.set_pane_color((0.,0.,0.,0.))
        def _axcolor(color):
            for spine in self.ax.get_children():
                if isinstance(spine, mpl.spines.Spine):
                    spine.set_color(color)
            self.ax.tick_params(color=color, labelcolor=color)
            self.ax.xaxis.label.set_color(color)
            self.ax.yaxis.label.set_color(color)
            self.ax.zaxis.label.set_color(color)
            mpl.rcParams['axes.labelcolor'] = color
        if self.guistyle == 'dark':
            self.fig.set_facecolor('black')
            self.ax.set_facecolor('black')
            _clearpanes()
            _axcolor('white')
            mpl.rcParams['text.color'] = 'white'
            mpl.rcParams['axes.prop_cycle'] = cycler(color=(
            'cyan', 'yellow', 'lime', 'red', 'azure', 'blue', 'gold',
            'yellowgreen', 'orange', 'blue', 'beige', 'fuchsia'))
        if not self.show_grid:
            self.ax.grid(False)
            _axcolor((0.,0.,0.,0.))
            _clearpanes()


    def _init(self, i):
        """
        Initialize the simulation by evaluating the physics engine for a specified number of steps.

        Args:
            i (int): The number of steps to initialize the simulation.

        Note:
            This method evaluates the physics engine for the specified number of steps to set up the simulation.
        """

        for _ in trange(i, desc='Evaluating motion for each frame', unit='frames'):
            self._engine.evaluate()
        tqdm.write('Calculations finished, Starting interactive window...')
    def _animate(self, i):
        """
        Animate the simulation by updating the visualization.

        Args:
            i (int): The frame index for animation.

        Note:
            This method updates the visualization of the simulation for a given frame index.
        """
        
        
        
        maxim = len(self._engine.bodies[0].pos.X.hist)
        ind = int(i*self._frameskip)
        step = self._plotskip
        while ind >= maxim:
            ind =- 1 
        self.ax.clear()
        self.ax.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')
        self.ax.set_box_aspect(None, zoom = self.zoom_slider.val)
        if not self.autoscale:
            self.ax.set_box_aspect((1,1,1), zoom = self.zoom_slider.val)
            self.ax.set_autoscale_on(False)
            if self.focus_body is not None:
                (limx, limy, limz) = (float(m) for m in self.focus_body.pos[ind])
                if self.focus_range is None:
                    self.focus_range = max((max(self.focus_body.pos - bod.pos) for bod in self._engine.bodies))
                if isinstance(self.focus_body.pos.X[ind], DecType) or isinstance(self.focus_range, DecType):
                    self.focus_range = float(self.focus_range)
                    limx, limy, limz = float(limx),float(limy), float(limz)
            else:
                limx, limy, limz = 0,0,0
                if self.focus_range is None:
                    self.focus_range = max(max(*bod.pos[ind]) for bod in self._engine.bodies)
            
           
            self.ax.set_xlim(xmin=(limx-self.focus_range),
                             xmax=(limx+self.focus_range))
            self.ax.set_ylim(ymin=(limy-self.focus_range),
                             ymax=(limy+self.focus_range))
            self.ax.set_zlim(zmin=(limz-self.focus_range),
                             zmax=(limz+self.focus_range))
        else:
            self.ax.set_autoscale_on(True)
            self.ax.set_autoscalez_on(True)
        for plane in self._engine.planes:
            xl, yl, zl, = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
            pl_const = np.array([[plane[1], plane[1]], [plane[1], plane[1]]])
            points = {'x':(pl_const,np.array([[yl[0],yl[1]],[yl[0],yl[1]]]),
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'y':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),pl_const,
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'z':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),
                           np.array([[yl[0],yl[0]],[yl[1],yl[1]]]),pl_const)}
            self.ax.plot_surface(*points[plane[0]], zorder=1,color=('xkcd:azure', 0.5), clip_on=False)
        for b in self._engine.bodies:
            _poshist = list(list(float(m) for m in _b.hist[0:ind:step]) for _b in (b.pos.X, b.pos.Y, b.pos.Z))
            if self.show_velocity or self.show_acceleration:
                _pos = [float(m) for m in b.pos[ind]]
                if self.show_velocity:
                    _vel = [float(m) for m in b.vel[ind]]
                    self.ax.quiver(*_pos, *_vel, length=self.vector_size, color='red',
                                   clip_on=False, zorder=8)
                if self.show_acceleration:
                    _acc = [float(m) for m in b.acc[ind]]
                    self.ax.quiver(*_pos, *_acc, length=self.vector_size, color='green',
                                   clip_on=False, zorder=8)
            self.ax.plot(*_poshist, label=f'{b.identity}', color=b.color,
                         zorder=7, clip_on=False)
            if self.body_model == 'dots':
                self.ax.scatter(*list(float(m) for m in b.pos[ind]),
                                marker='o', zorder=4, clip_on=False, color=b.color)
            if self.body_model == 'wireframe':
                self.ax.plot_wireframe(*sphere(b.pos[ind], b.radius.c()), clip_on=False,
                                       zorder=2, color=b.color)
            if self.body_model == 'surface':
                self.ax.plot_surface(*sphere(b.pos[ind], b.radius.c()), clip_on=False,
                                     zorder=2, color=b.color)
            if self.labelling_type == 'label':
                self.ax.text(*b.pos[ind], b.identity, zorder=10)
            if self.show_shadows:
                self.ax.plot(*_poshist[0:2],[(self.focus_body.pos.Z[ind]-self.focus_range)]*
                             len(_poshist[2]),
                    color='black', zorder=1.5, clip_on=False)
        if self.labelling_type == 'legend':
            self.ax.legend()


    def start(self, eval_length=None, fps=None, frameskip=1, plotskip=1):
        """
        Start the animation of the simulation, and initialise output window.

        Args:
            frames (int, optional): The total number of frames for the animation.
            interval (int, optional): The interval between frames in seconds.
            duration (int, optional): The total duration of the animation in seconds.
            fps (int, optional): The frames per second for the animation.

        Note:
            This method starts the animation of the simulation based on the provided animation parameters.
        """
        self._frameskip, self._plotskip = frameskip, plotskip
        f,inv = eval_length, (1/fps)/1000

        print('Starting Simulation Instance, Running Calculations:')
        anim = animation.FuncAnimation(self.fig, func = self._animate,
                                       init_func = self._init(f), interval=inv, frames=f) 
        plt.show()
#END of Simulation Class       
class SolarSystemMB(Simulation):
        def __init__(self, dt=1000, show_grid: bool = True,
                     show_shadows: bool = False,
                     show_acceleration: bool = False,
                     show_velocity: bool = False,
                     vector_size: int | float = 1,
                     labelling_type: str = 'label',
                     body_model: str = 'dots',
                     guistyle: str = 'dark'):
            name = 'Major Bodies in Solar System'
            engine = PhysEngine(dt)
            bodies = list(horizons_query(obj_id) for obj_id in (
                '10', '199', '299','399', '499', '599', '699', '799', '899'))
            engine.attach_bodies(bodies)
            engine.make_relative_to(bodies[0])
            focus_body, focus_range, autoscale = bodies[0], None, False 
            super().__init__(name, engine, focus_body,
                             focus_range, autoscale, show_grid,
                             show_shadows, show_acceleration, show_velocity,
                             vector_size, labelling_type, body_model, guistyle)
            #print('Note: the suggested start() Parameters are frames=10000<->50000, fps=30, frameskip=100, plotskip=50.')
        def start(self, eval_length=50000, fps=30, frameskip=1000, plotskip=500):
            super().start(eval_length, fps, frameskip, plotskip)