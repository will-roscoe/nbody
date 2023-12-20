#### 3rd Party Libs/Packages
import numpy as np
# ⤤ for saving engine as binary
from tqdm import tqdm, trange
# ⤤ for progress bars so then user has expected time to completion per part
#### Local imports
from ..tools import errors as e
# ⤤ standardized error messages
from .base import typecheck, _O, Iterable, NumType, Vector, VectorType, NullVector
from .body import Body, G
# ⤤ neccesary imports for engine



class Engine:
    '''This object acts as the physics engine, and simulates the gravitational effects of bodies, and can also\
        simulate collisions as well.
### Parameters
|Parameter| Required |Type| Description|
|---|---|---| ---|
|`dt` | ✕ | `base.NumType` | interval in seconds between each datapoint in the simulation.\
    Lowest will produce most accurate results. Default is `1`(second). |
|`checking_range` | ✕ | `int` | range of points to check for collisions over for each point. Default is `3`.|

### Attributes
|Attribute| Type| Description|
|---|---|---|
| `self.bodies` |`list`  | internal collection of bodies that have been attached. |
| `self.planes` |`list`  | internal collection of planes that have been defined. |
| `self.fields` |`list`  | internal collection of fields that have been defined. |
| `self.dt` |`base.NumType`  | see `dt` above |
| `self.do_collisions` |`bool`  | whether to calculate the effect of collisions with planes or other bodies when\
    simulating the trajectories. Default=`True`. |
| `self.do_bodygravity` |`bool`  | whether to calculate the effect of gravity due to other bodies when simulating the\
    trajectories. Default=`True`.|
| `self.do_fieldgravity` |`bool`  | whether to calculate the effect of fields when simulating the trajectories.\
    Default=`True`.|

### Usage
#### `attach_bodies(new_bodies)`
   - attaches `core.Body` objects to the `Engine` instance. `new_bodies` must be passed as a `base.Iterable`.

#### `save_as(dump='engine',file_name='nbody_data')`
   - saves the current state and all child objects as an .npz file.
        * `dump` = either  `engine` or `bodies`, choosing to save engine with all children or just the\
              `core.Bodies` objects.
#### `load_as(objects='engine',file_name='nbody_data')`
   - loads the chosen file data into the `Engine` instance. `objects` = either `engine` or `bodies` and this specifies\
      the target object(s) to load.
#### `make_relative_to(target_body)`
   - reinitialises all bodies in the instance, changing the initial conditions so that `target_body` is located at\
    (0,0,0) and the rest of the objects have the same initial conditions relative to `target_body`.
#### `create_acceleration(accel_vector)`
   - creates a constant field of acceleration that can affect the simulated trajectories of bodies.
     * `accel_vector` must be an iterable or `base.VectorType` instance.
>[!WARNING]
> if the value of an acceleration field or a body's velocity is sufficiently high, a body may fail to register\
      a collision with a plane or another body.
>
> In these instances, it is neccesary to reduce `dt`. *(and increase `checking_range`, which will sample more \
    points resulting in an earlier check collision being registered.)*
#### `create_plane(const_axis= 'z', const_val= 0)`
   - creates an infinite plane with 0 thickness, parallel to a chosen euclidean plane. `const_axis` specifies\
    the plane orientation, e.g, `'z'` would return a plane parallel to the x-y plane. `const_val` specifies\
        the coordinate value it should be locked to.
#### `simulate(self,intervals)`
  - runs a simulation over a specified amount of `intervals` (must be `int`). does not return anything,\
    but each of the bodies will contain the new data, and the `Engine` instance can now be passed on to a \
        `Visual` object to create an output.

    '''
    def __init__(self,dt=1,checking_range=3):
        self.bodies, self.planes, self.fields= [], [], []
        (self.dt,self._rangechk) = typecheck(((dt,NumType),(checking_range,NumType)))
        self.do_collisions,self.do_bodygravity,self.do_fieldgravity = True,True,True
    
    def __len__(self):
        if len(self.bodies) is not None:
            # returning length of data for body 0, assuming all are the same.
            return len(self.bodies[0].pos)
        else:
            return 0
    
    def attach_bodies(self, new_bodies):
        '''attaches `core.Body` objects to the `Engine` instance. `new_bodies` must be passed as a `base.Iterable`.
        '''
        if isinstance(new_bodies,Iterable):
            for new_body in new_bodies:
                typecheck((new_body, Body))
                self.bodies.append(new_body)
            tqdm.write(f'«Engine» → {len(self.bodies)} bodies currently attached to engine.')
    
    
    def _loadeng(self,eng):
        # used for loading an engine from an npz file.
        self = eng #noqa: F841
    
    
    def save_as(self,dump='engine',file_name='nbody_data'):
        '''saves the current state and all child objects as an .npz file.
        * `dump` = either  `engine` or `bodies`, choosing to save engine
        with all children or just the `core.Bodies` objects.'''
        _saveobjs = {'bodies':{'bodies':self.bodies},'engine':{'engine':self}}
        # saving object as npz file
        np.savez(f'{file_name}.npz', **_saveobjs[dump])

    
    def load_as(self,objects='engine',file_name='nbody_data'):
        '''loads the chosen file data into the `Engine` instance. `objects` = either
        `engine` or `bodies` and this specifies the target object(s) to load.'''  
        _loadobjs = {'bodies':("self.attach_bodies(objs['bodies'])",),
                     'engine':("self._loadeng(objs['engine'])",)}
        
        objs = np.load(f'{file_name}.npz',allow_pickle=True) #noqa
        if objects == 'engine':
            self = objs['engine']
        elif objects == 'bodies':
            self.bodies = objs['bodies']
        else:
            raise ValueError
        objs.close()
    
    def make_relative_to(self,target_body):
        '''reinitialises all bodies in the instance, changing the initial conditions so that `target_body` is located
        at (0,0,0) and the rest of the objects have the same initial conditions relative to `target_body`.'''
        for body in self.bodies:
            if body != target_body:
                body._reinitialise((body.pos-target_body.pos), (body.vel-target_body.vel))
                body.identity = f'{body.identity}:Rel.{target_body.identity}'
        
        target_body._reinitialise((0,0,0),(0,0,0))
        tqdm.write(f"«Engine» → Bodies' positions and velocities have been made relative to {target_body.identity}.")
        target_body.identity = f'{target_body.identity}(Static)'
  

    
    def create_acceleration(self, accel_vector):
        '''creates a constant field of acceleration that can affect the simulated trajectories of bodies.
     * `accel_vector` must be an iterable or `base.VectorType` instance.'''
        _acc = _O(accel_vector)
        if len(_acc) == 3 and isinstance(_acc[0], NumType):
            # adding a vector instance to the fields list
            self.fields+= [Vector(li=_acc)] 
            tqdm.write(f'«Engine» → constant acceleration {accel_vector} has been initialized.')
        else:
            e.raise_type_error('accel_vector', (Iterable, *VectorType), accel_vector)
    
    def create_plane(self, const_axis ='z', const_val  = 0):
        '''creates an infinite plane with 0 thickness, parallel to a chosen euclidean plane. `const_axis` specifies the
        plane orientation, e.g, `'z'` would return a plane parallel to the x-y plane. `const_val`
        specifies the coordinate value it should be locked to.'''
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, NumType):
            self.planes.append([const_axis, const_val])
            tqdm.write(f'«Engine» → constant plane {const_axis}={const_val} has been initialized.')
        else:
            e.raise_value_error('const_axis,const_val',(str,*NumType),(const_axis,const_val))

       
    def _check_collision(self,body,co_restitution=0):
        if self.do_collisions == True:
            # has it found a collision?
            returned_coll = False
            for body2 in self.bodies:
                if body != body2: 
                    # get distance between bodies
                    body_dist = (body2.pos - body.pos).magnitude()             
                    # check if distance is less than radii
                    if body_dist < (body2.radius + body.radius).c():
                        returned_coll = True
                        n = (body2.pos-body.pos).unit()
                        meff = 1/((1/body2.mass)+(1/body.mass))
                        rel_v = body2.vel - body.vel
                        impulse =  n *(n*rel_v)*(1 + co_restitution) * meff
                        dv = impulse/body.mass
                        # calc vel change: dv = n^ * [n^*[v2 -v1]*[1+e]*m"]/m1, n^ is normal unit vector,
                        # m" is effective mass as above.
                        #dv = n * -1*(( (  n   *(body.vel-body2.vel) )*(   1+co_restitution)*meff     )/body.mass.c())
                        return (dv+body.vel), False
            
            unitcomps = {'x':Vector((1,0,0)),'y':Vector((0,1,0)),'z':Vector((0,0,1))}
            for pl in self.planes:
                # check for plane collision. we estimate a range of times as the plane has no thickness,
                # and time interval would have to be extremely miniscule compared to velocity to always catch.
                body_dists = list(abs(((body.pos + body.vel*m*self.dt)[pl[0]] - pl[1])) for m in range(self._rangechk))
                if any([(body_dists[i] < body.radius) for i in range(self._rangechk)]):
                    on_plane = ((body_dists[0]/body.radius <= 1.01) if len(body_dists) == 1 
                                else (body_dists[0]/body.radius <= 1.01 and body_dists[-1]/body.radius <= 1.01))
                    returned_coll = True
                    # similar to above, but normals are precalculated.
                    return (unitcomps[pl[0]]*(-2*(body.vel*unitcomps[pl[0]])) + body.vel)*co_restitution, on_plane
            
            if not returned_coll:
                return 1*body.vel, False
        else:
            return 1*body.vel, False
    
    
    def _find_gravity(self,body):
        res = NullVector()
        if self.do_bodygravity == True:
            for bod2 in self.bodies: 
                if bod2 != body:
                    dist_12 = body.pos - bod2.pos
                    unit_12,mag_12 = dist_12.unit(),dist_12.magnitude()
                    force_on_bod1_t = (-1*G * body.mass * bod2.mass)
                        # f = -G * M1 * M2
                    force_on_bod1 = unit_12*(force_on_bod1_t/(mag_12)**2)
                        # F = r^ * f/(|r|**2)
                    # sum all forces
                    res += force_on_bod1
        
        return res/body.mass
    
    
        
    
    def simulate(self,intervals):
        '''runs a simulation over a specified amount of `intervals` (must be `int`). does not return anything, but each
          of the bodies will contain the new data, and the `Engine` instance can now be passed 
          on to a `Visual` object to create an output.'''
        if isinstance(intervals, int) and len(self.bodies) > 0:
            for _ in trange(intervals, 
                            desc=f'«Engine» → Evaluating motion for each interval of {self.dt} seconds', unit='ints'):
                _temp = [[0,0,0] for _ in self.bodies]
                for i,body in enumerate(self.bodies):
                    _temp[i] = [*self._check_collision(body, body.bounce),self._find_gravity(body)]
                
                for i,body in enumerate(self.bodies):
                    col_vel, on_plane, acc_g = _temp[i]
                    if not on_plane and self.do_fieldgravity:
                        fieldvel = list(sum(v.c(i) for v in self.fields) for i in range(3))
                    else:
                        fieldvel = NullVector()
                    body.update(dt=self.dt,vel_next=(col_vel+fieldvel).c(),acc_change=acc_g)
            
            if intervals+1 == len(self.bodies[0].pos):
                tqdm.write(
        f'«Engine» → Finished Evaluating {intervals} intervals, ~{len(self.bodies[0].pos)} total intervals.',
                )
    
    def barycenter(self, index):
        '''gets the position of the center of mass of the engine'''
        mass_dist = NullVector()
        for b in self.bodies:
            mass_dist += Vector(b.pos[index])*b.mass # add up all mass*pos
        total_mass =sum([b.mass.c() for b in self.bodies]) 
        return mass_dist/total_mass # divide by mass
    
    def __getitem__(self, ind):
        if ind == '_data_':
            return {
                'dt':self.dt,
                'chr':self._rangechk,
                'docoll':self.do_collisions,
                'dobg':self.do_bodygravity,
                'dofg':self.do_fieldgravity,
                'planes':self.planes,
                'fields':self.fields,
                'bodies':self.bodies,
            }
            
            