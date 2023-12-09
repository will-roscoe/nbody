import numpy as np
from tqdm import tqdm, trange
from decimal import Decimal

from ..tools import errors as e
from .base import typecheck, _O, Iterable, NumType, Vector, VectorType, NullVector, DecType
from .body import Body, G



class Engine:
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
        if isinstance(new_bodies,Iterable):
            for new_body in new_bodies:
                typecheck((new_body, Body))
                self.bodies.append(new_body)
            tqdm.write(f'«Engine» → {len(self.bodies)} bodies currently attached to engine.')
    
    
    def _loadeng(self,eng):
        # used for loading an engine from an npz file.
        self = eng
    
    
    def save_as(self,dump='engine',file_name='nbody_data'):
        _saveobjs = {'bodies':{'bodies':self.bodies},'engine':{'engine':self}}
        # saving object as npz file
        with open(f'{file_name}.npz','wb') as file:
            np.savez(file, **_saveobjs[dump])


    def load_as(self,objects='engine',file_name='nbody_data'):   
        _loadobjs = {'bodies':("self.attach_bodies(objs['bodies'])",),
                     'engine':("self._loadeng(objs['engine'])",)}
        
        with open(f'{file_name}.npz','rb') as file:
            objs = np.load(file,allow_pickle=True)
            for func in _loadobjs[objects]:
                try:
                    # usage limited to values in dict above.
                    eval(func)
                except KeyError:
                    raise LookupError(f'cannot find {objects} value in "{file}"')


    def make_relative_to(self,target_body):
        for body in self.bodies:
            if body != target_body:
                body._reinitialise((body.pos-target_body.pos), (body.vel-target_body.vel))
                body.identity = f'{body.identity}:Rel.{target_body.identity}'
        
        target_body._reinitialise((0,0,0),(0,0,0))
        tqdm.write(f"«Engine» → Bodies' positions and velocities have been made relative to {target_body.identity}.")
        target_body.identity = f'{target_body.identity}(Static)'
  

    
    def create_acceleration(self, accel_vector):
        _acc = _O(accel_vector)
        if len(_acc) == 3 and isinstance(_acc[0], NumType):
            # adding a vector instance to the fields list
            self.fields+= [Vector(li=_acc)] 
            tqdm.write(f'«Engine» → constant acceleration {accel_vector} has been initialized.')
        else:
            e.raise_type_error('accel_vector', (*Iterable, *VectorType), accel_vector)
    
    def create_plane(self, const_axis ='z', const_val  = 0):
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, NumType):
            self.planes.append([const_axis, const_val])
            tqdm.write(f'«Engine» → constant plane {const_axis}={const_val} has been initialized.')
        else:
            e.raise_value_error('const_axis,const_val',(str,*NumType),(const_axis,const_val))

       
    def _check_collision(self,body,co_restitution=0):
        if self.do_collisions == True:
            # has it found a collision?
            returned_coll = False
            for bod in self.bodies:
                if body != bod: 
                    # get distance between bodies
                    body_dist = (Vector(bod.pos.c()) - Vector(body.pos.c())).magnitude()             
                    # check if distance is less than radii
                    if body_dist < (bod.radius.c() + body.radius.c()):
                        (returned_coll,n,meff) = (True,Vector(bod.pos-body.pos)/body_dist,1/((1/bod.mass.c())+(1/body.mass.c())))
                        #calc vel change: dv = n^ * [n^*[v2 -v1]*[1+e]*m"]/m1, n^ is normal unit vector, m" is effective mass as above.
                        dv = n*(-((n*(body.vel - bod.vel))*(1+co_restitution)*meff)/body.mass.c())
                        return (dv+body.vel), False
            
            unitcomps = {'x':Vector((1,0,0)),'y':Vector((0,1,0)),'z':Vector((0,0,1))}
            for pl in self.planes:
                # check for plane collision. we estimate a range of times as the plane has no thickness, and time interval would have to be extremely miniscule compared to velocity to always catch.
                body_dists = list(abs((Vector(body.pos.c()) +
                Vector(body.vel.c())*m*self.dt)[pl.c] - pl.v) for m in range(self._rangechk))
                if any([(body_dists[i] < body.radius.c()) for i in range(self._rangechk)]) or\
                body_dists[0] < body.radius.c():
                    on_plane = (True if body_dists[0]/body.radius.c() <= 1.01 and body_dists[-1]/body.radius.c() <= 1.01 else False)
                    returned_coll = True
                    # similar to above, but normals are precalculated.
                    return (unitcomps[pl.c]*(-2*(body.vel*unitcomps[pl.c])) + body.vel)*co_restitution, on_plane
            
            if not returned_coll:
                return Vector(body.vel.c()), False
        else:
            return Vector(body.vel.c()), False
    
    
    def _find_gravity(self,body):
        res = NullVector()
        if self.do_bodygravity == True:
            for bod2 in self.bodies: 
                if bod2 != body:
                    dist_12 = body.pos - bod2.pos
                    unit_12,mag_12 = dist_12.unit(),dist_12.magnitude()
                    
                    if isinstance(body.mass.c(),DecType) or isinstance(bod2.mass.c(),DecType):
                        b1m, b2m = Decimal(body.mass.c()), Decimal(bod2.mass.c())
                        force_on_bod1_t = (-1*Decimal(G) * b1m * b2m)
                        force_on_bod1 = unit_12*(force_on_bod1_t/Decimal(mag_12)**2)
                    else:
                        force_on_bod1_t = (-1*G * body.mass.c() * bod2.mass.c())
                        # f = -G * M1 * M2
                        force_on_bod1 = unit_12*(force_on_bod1_t/(mag_12)**2)
                        # F = r^ * f/(|r|**2)
                    # sum all forces
                    res += force_on_bod1
        
        return res/body.mass.c()

    
    def simulate(self,intervals):
        if isinstance(intervals, int) and len(self.bodies) > 0:
            for _ in trange(intervals, desc=f'«Engine» → Evaluating motion for each interval of {self.dt} seconds', unit='ints'):
                _temp = [[0,0,0] for _ in self.bodies]
                
                for i,body in enumerate(self.bodies):
                    _temp[i] = [*self._check_collision(body, body.bounce),self._find_gravity(body)]
                
                for i,body in enumerate(self.bodies):
                    col_vel,on_plane,acc_g = _temp[i]
                    if not on_plane and self.do_fieldgravity:
                        fieldvel = list(sum(v.c(i) for v in self.fields) for i in range(3))
                    else:
                        fieldvel = NullVector()
                    body.update(dt=self.dt,vel_next=(col_vel+fieldvel).c(),acc_change=acc_g)
            
            if intervals+1 == len(self.bodies[0].pos):
                tqdm.write(f'«Engine» → Finished Evaluating {intervals} intervals, ~{len(self.bodies[0].pos)} total intervals.')

    def barycenter(self, index):
        mass_dist = NullVector()
        for b in self.bodies:
            mass_dist += Vector(b.pos[index])*b.mass.c()
        total_mass =sum([b.mass.c() for b in self.bodies]) 
        return mass_dist/total_mass
    
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
                'bodies':self.bodies
            }
            
            