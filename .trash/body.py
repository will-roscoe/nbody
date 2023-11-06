from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta
import astropy.units as u
import re

from core import HistoricVariable, HistoricVector, Vector
import errors as e


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
        if isinstance(mass, (float, int)):
            self.mass = HistoricVariable(mass,
                                identity=f'{self.identity}_mass',
                                units='kg')
        else:
            e.raise_type_error('mass', (int, float), mass)
        if isinstance(radius, (float, int)):
            self.radius = HistoricVariable(radius,
                                identity=f'{self.identity}_rad',
                                units='m')
        else:
            e.raise_type_error('radius', (int, float), radius)


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
        This method updates the position, velocity, and acceleration of a body over a small time interval `dt`.
        It ensures that the lengths of position and velocity vectors are consistent and calculates the acceleration.

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
                acc_change: list | tuple=None, ) -> None:
        """
        Update the position and velocity of a body over a small time interval.

        This method updates the position and velocity of a body based on the time step `dt`, velocity changes (`vel_change`),
        and acceleration changes (`acc_change`). It ensures that the lengths of position and velocity vectors are consistent.

        Args:
            dt (int or float, optional): The time step for the update (default is 1).
            vel_change (list or tuple, optional): Velocity changes as a 3-element list or tuple (default is None).
            acc_change (list or tuple, optional): Acceleration changes as a 3-element list or tuple (default is None).

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
        if len(self.pos) - len(self.vel) != 0:
            self.evaluate(dt)
        if len(self.pos) - len(self.vel) == 0:
            if acc_change is not None:
                if isinstance(acc_change, (list, tuple)):
                    if len(acc_change) == 3:
                        self.vel.next(self.vel + Vector((self.acc + acc_change))*dt)
                    else:
                        e.raise_component_error('acc_change', acc_change)
                else:
                    e.raise_type_error('acc_change', (list, tuple), acc_change)
            else:
                self.vel.next(self.vel + Vector((self.acc*dt)))
            if vel_change is not None:
                if isinstance(vel_change, (list, tuple)):
                    if len(vel_change) == 3:
                        self.pos.next(self.pos + Vector((self.vel + vel_change))*dt)
                    else:
                        raise e.raise_component_error('vel_change', vel_change)
                else:
                    e.raise_type_error('vel_change', (list, tuple), vel_change)
            else:
                self.pos.next(self.pos + Vector((self.vel*dt)))
        else:
            e.raise_evaluation_error((self.pos, self.vel))
#END of Body Class



def make_horizons_object(searchquery, observer='0', time='2023-11-03'):
    """
    Create a Horizons object representing a celestial body using Horizons query data.

    This function queries NASA's Horizons system for information about a celestial body and creates a
    `Body` object with the retrieved data, including mass, initial position, initial velocity, radius,
    and identity.

    Args:
        searchquery (str): The identifier or name of the celestial body to query.
        observer (str, optional): The observer location (default is '0' for the solar system barycenter).
        time (str, optional): The date for which to retrieve data in the 'YYYY-MM-DD' format (default is '2023-11-03').

    Returns:
        Body: A `Body` object representing the celestial body with mass, initial position, initial velocity,
        radius, and identity information.

    Raises:
        LookupError: If mass or radius information could not be found in the query output.
    """
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
    _object = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors(get_raw_response=True)
    _object_tab = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors()
    _m = re.search(r"\s+?Mass\s+?([\S]+).*?=\s+?([\S]+?)\s+?", _object)
    _r = re.search(r'\s+?radius\s+?([\S]+).*?=\s+?([\S]+?)\s+?', _object)
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _object)
    if _m:  
        _mbase, _mexp_s = _m.groups()[1], _m.groups()[0]
        mass = float(_mbase.split('+-')[0]) * 10**(float(_mexp_s.split('x10^')[1]))
    else: 
        raise LookupError('could not find Mass in output')
    if _r:
        _rad =_r.groups()[1]
        radius = float(_rad.split('+-')[0])*1000
    else:
        raise LookupError('could not find radius (Vol. mean Radius) in output')
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        print('could not find name for object, reverting to placeholder name.')
        name = None
    x, y, z = [_object_tab[pos].quantity[0].to_value(u.m) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [_object_tab[pos].quantity[0].to_value(u.m/u.s) for pos in ('vx', 'vy', 'vz')]
   
    return Body(mass=mass, init_pos=(x,y,z), init_vel=(vx,vy,vz), radius=radius, identity=name)