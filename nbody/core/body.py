import math
# ⤤ sqrt and pi
from ..tools import errors as e
# ⤤ standard error messages
from .base import (typecheck, _O, _V, 
                   Iterable, NumType, NoneType, VectorType,
                   HistoricVector, Variable, Vector, NullVector)

try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)
# ⤤ if scipy isnt found, revert to an approximation


class Body:
    '''This object is the representation of a single (solid) body with finite mass and size.
Should be placed into an `core.Engine` instance to perform a simulation.

### Parameters
|Parameter| Required |Type| Description|
|---|---|---| ---|
|`mass` | ✓ | `base.NumType` | mass of the body in kg. |
|`init_pos` | ✓ | `base.Iterable` or `base.VectorType` | initial position of the body in m.|
|`init_vel` | ✕ | `base.Iterable` or `base.VectorType` | initial velocity of the body in ms^-1. Default is `(0,0,0)`.|
|`radius` | ✕ | `base.NumType` | radius of the body in m, modelled as a sphere. Default is `0`. if `radius=0`,\
collision physics will not be simulated, and will act as a point particle. |
|`bounce` | ✕ | `base.NumType` | coefficient of restitution for the body. the factor describing how much energy\
is lost in a collision. `1` will give perfect elastic collision, `0` will give a perfect inelastic collision.\
    Default is `0.999`.  |
|`color` | ✕ | `str` | color of the body when visualised. can used colors defined in matplotlib docs.|
### Attributes
|Attribute| Type| Description|
|---|---|---|
| `self.identity` |`str`  | identifiable string to represent the body by. |
| `self.mass` |`base.Variable`  | mass of the body. |
| `self.pos` |`base.HistoricVector`  | positional vectors of body. |
| `self.vel` |`base.HistoricVector`  | velocity vectors of body. |
| `self.acc` |`base.HistoricVector`  | acceleration vectors of body. |
| `self.radius` |`base.Variable`  | radius of the body. |
| `self.bounce` |`base.NumType`  | coefficient of restitution used when the body is in a collision. |
| `self.color` |`str`  | color of the body when visualised. |
### Usage
#### `update(dt=1,vel_change=None,acc_change=None,vel_next=None)`
   - evaluates the bodies current position, velocity and acceleration, along with any optional changes,
   \over a time interval `dt`.

#### `_reinitialise(init_pos=None,init_vel=None)`
  - resets the body to an initial position and length, erasing all other datapoints for position, 
  velocity and acceleration.

#### `get_(item, ind=-1, plotskip=0, c_mass=None, engine=None)`
   - calculates and returns a physical quantity of the `Body` object. either `period`, `sma`(semi major axis), or 
   `ke` (kinetic energy). 
### Indexing
- Supports numerical indexing but not slicing.

| Index | Returns|
|---|---|
|`'pos'`| all points of position as a 2d array.|
|`'vel'`| all points of velocity as a 2d array.|
|`'acc'`| all points of acceleration as a 2d array.|
|`'current'`| 2d array containing the current position, velocity and acceleration.|
|`'info'`| dictionary containing all non kinematic properties of the body|
|`'x'`| 2d array containing all x components of position, velocity and acceleration vectors.|
|`'y'`| 2d array containing all y components of position, velocity and acceleration vectors.|
|`'z'`| 2d array containing all z components of position, velocity and acceleration vectors.|
    '''
    def __init__(self,mass,init_pos,init_vel=(0,0,0),
                radius=0,bounce=0.999,color=None,identity=None):

        if isinstance(identity,str):
            self.identity = identity
        elif identity is None:
            self.identity = f'body[m={mass}]'
        else:
            e.raise_type_error('identity', str, identity)
        
        init_pos, init_vel = _O(init_pos), _O(init_vel)
        typecheck(((init_pos,Iterable),(init_vel,Iterable),(mass,NumType),
                     (radius,NumType),(bounce,NumType),(color,(NoneType, str))))
        
        _id = self.identity
        
        self.pos = HistoricVector(li=init_pos,identity=f'{_id}_pos',units_v='metre')
        self.vel = HistoricVector(li=init_vel,identity=f'{_id}_vel',units_v='metre / second')
        self.acc = HistoricVector(0,0,0,identity=f'{_id}_acc',units_v='metre / second^2')
        
        self.mass = Variable(mass,identity=f'{_id}_mass',units='kg')
        self.radius = Variable(radius,identity=f'{_id}_rad',units='metre')
        self.bounce = bounce
        self.color = color
    
    
    def __str__(self):
        return f'Body("{self.identity}",\n\
    mass="{self.mass.c()} {self.mass.units}",\n\
    currentpos="{self.pos.c()} {self.pos.units}",\n\
    currentvel="{self.vel.c()} {self.vel.units}",\n\
    currentacc="{self.acc.c()} {self.acc.units}")'
    
    
    def __repr__(self):
        return f'Body("{self.identity}", m={self.mass.c()}, r={self.pos.c()},\
v={self.vel.c()}), a={self.acc.c()})'
    
    def __getitem__(self, ind):
        if isinstance(ind, int):
            return (self.pos[ind], self.vel[ind], self.acc[ind])
        elif isinstance(ind, str):
            return {'pos':self.pos['all'], 'vel': self.vel['all'], 'acc':self.acc['all'], 
    'current':list(d.c() for d in (self.pos, self.vel, self.acc)),'info':{'identity':self.identity, 'mass':self.mass, 
    'radius':self.radius, 'color':self.color, 'bounce':self.bounce},'x':list(d.X for d in 
                                                                             (self.pos, self.vel, self.acc)),
    'y':list(d.Y for d in (self.pos, self.vel, self.acc)),'z':list(d.Z for d in (self.pos, self.vel, self.acc)),
    '_data_':{'pos':self.pos['all'], 'vel': self.vel['all'], 'acc':self.acc['all'],'identity':self.identity, 
    'mass':self.mass, 'radius':self.radius, 'color':self.color, 'bounce':self.bounce}}[ind]
    
    def update(self,dt=1,vel_change=None,acc_change=None,vel_next=None):
        '''evaluates the bodies current position, velocity and acceleration, along with any optional changes,
over a time interval `dt`. 
    `dt:Number` - interval to calculate over.
    `vel_change,acc_change,vel_next: Iterable | VectorType` - changes to the bodies current condition.
        '''
        vel = _V(vel_next) if vel_next else self.vel
        
        if acc_change is not None:
            acc_change = _V(acc_change)
            self.vel.next((acc_change*dt + vel))
            self.acc.next(acc_change.c())

        else:
            self.vel.next((self.acc*dt)+vel)
            self.acc.next(NullVector())
    
        if vel_change is not None:
            vel_change = _V(vel_change)
            self.pos.next(self.pos + (vel + vel_change)*dt)
        else:
            self.pos.next(self.pos + vel*dt)

    
    
    def _reinitialise(self,init_pos=None,init_vel=None):
        '''resets the body to an initial position and length, erasing all other datapoints for position,
        velocity and acceleration.
        '''
        self.acc = HistoricVector(0,0,0,identity=f'{self.identity}_acc',units_v=self.acc.units)
        
        if init_pos != None:
            if isinstance(init_pos,(Iterable,*VectorType)):
                self.pos = HistoricVector(li=init_pos,identity=f'{self.identity}_pos',units_v=self.pos.units)
            else:
                e.raise_type_error('init_pos',Iterable,init_pos)
        
        if init_vel != None:
            if isinstance(init_vel,(Iterable,*VectorType)):
                self.vel = HistoricVector(li=init_vel,identity=f'{self.identity}_vel',units_v=self.vel.units)
            else:
                e.raise_type_error('init_vel',(Iterable,*VectorType),init_vel)

    
    
    def get_(self, item, ind=-1, plotskip=0, c_mass=None, engine=None):
            '''calculates and returns a physical quantity of the `Body` object. 
            either `period`, `sma`(semi major axis), or `ke` (kinetic energy).
            '''
            def sma():
                if not plotskip >= ind:
                    a = max((self.pos[i]-engine.barycenter(ind)).magnitude() for i in range(0,ind,plotskip))
                    return a
                else:
                    return 'NaN'
            def per():
                a = sma()
                if a != 'NaN':
                    return 2*math.pi*math.sqrt((a**3)/(G*(c_mass+self.mass).c()))
                else:
                    return a
            def ke():
                return ((Vector(self.vel[ind]).magnitude()**2)*(1/2)*self.mass).c()
            
            _get_lookup = {**dict.fromkeys(['sma', 'semi_major_axis'],sma),
                       **dict.fromkeys(['per', 'period'],per),
                       **dict.fromkeys(['ke', 'kinetic_energy'], ke)}
            return _get_lookup[item]()
        