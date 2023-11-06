import math


import errors as e



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

NoneType = type(None)

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
        elif all(isinstance(var, (int, float)) for var in (x, y, z)):
            self.X, self.Y, self.Z = x, y, z
        else:
            e.raise_list_type_error('l,x,y,z', (tuple, list, int, float), (li,x,y,z))

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
        elif (isinstance(x, (int, float)) and
              isinstance(y, (int, float)) and
              isinstance(z, (int, float))):
            self.X = HistoricVariable(x, f'{self.identity}_x', self.units)
            self.Y = HistoricVariable(y, f'{self.identity}_y', self.units)
            self.Z = HistoricVariable(z, f'{self.identity}_z', self.units)
        else:
            e.raise_type_error('l,x,y,z', (tuple, list, float, int), (li,x,y,z))

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
        return math.sqrt(sum([n**2 for n in self.c()]))
    def unit(self) -> tuple:
        """
        Get the unit vector in the same direction as the historic vector.

        Returns:
            tuple: A tuple representing the unit vector.
        """
        if float(self.magnitude()) == 0.:
            return (0,0,0)
        else: 
            return list(n/self.magnitude() for n in self.c())
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
