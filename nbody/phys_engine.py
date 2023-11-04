
try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)



from body import Body
from core import Vector, HistoricVector
import errors as e

#START of PhysEngine class
class PhysEngine:
    def __init__(self, dt: int | float = 1):
        self.bodies = []
        if isinstance(dt, (int, float)):
            self.dt = dt
        self.planes = []
    def attach_bodies(self, new_bodies):
        if isinstance(new_bodies, (list, tuple)):
            for i, new_body in enumerate(new_bodies):
                if isinstance(new_body, Body) or issubclass(type(new_body), Body):
                    self.bodies.append(new_body)
                else:
                    e.raise_type_error(f'new_body at index {i}', type(Body), new_body)
        else:
            e.raise_type_error('new_bodies', (list, tuple), new_bodies)
        print(f'{len(self.bodies)} bodies attached')

    def create_plane(self, const_axis='z', const_val = 0):
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, (int, float)):
            self.planes.append([const_axis, const_val])
        else:
            e.raise_value_error('const_axis,const_val',((str),(int,float)),(const_axis,const_val))
    def check_collision(self,body):
        unitcomps = {'x':Vector((1,0,0)), 'y':Vector((0,1,0)), 'z':Vector((0,0,1))}
        for bod in self.bodies:
            if body != bod: 
                temp_dist = body.pos - bod.pos
                if (Vector(temp_dist).magnitude() <= bod.radius.c() or
                Vector(temp_dist).magnitude() <= body.radius.c()):
                    del_v_bods=( Vector( Vector( body.vel - bod.vel ) /
                    ( ( 1/ body.mass.c() ) + ( 1/ bod.mass.c() ))) * -2)
                    break
        else:
            del_v_bods = Vector((0,0,0))
        for [pl_ax, pl_val] in self.planes:
            pl_norm = unitcomps[pl_ax]
            if body.pos[pl_ax] + pl_val <= body.radius:
                del_v_pls = body.vel + pl_norm*(-2*(body.vel*pl_norm))
                break
        else:
            del_v_pls = Vector((0,0,0))
        return (del_v_pls + del_v_bods).c()

    def find_gravity(self):
        temp = [[i, bod, HistoricVector(0,0,0)] for i, bod in enumerate(self.bodies)]
        for bod1_data in temp: 
            ind1, bod1, f1 = bod1_data[0], bod1_data[1], bod1_data[2]
            for bod2_data in temp: 
                ind2, bod2= bod2_data[0], bod2_data[1]
                if ind1 != ind2:
                    temp_dist = bod1.pos - bod2.pos
                    dist_12 = Vector(temp_dist)
                    unit_12 = dist_12.unit()
                    force_on_bod1_t = (-1*G * bod1.mass.c() * bod2.mass.c())
                    force_on_bod1 = unit_12*(force_on_bod1_t/(dist_12.magnitude())**2)
                    f1 = HistoricVector(li=(f1 + force_on_bod1))
                    bod1_data[2] = f1
        for bod_data in temp:
            bod, f = bod_data[1], bod_data[2]
            bod.acc.next(f/bod.mass.c())

    def evaluate(self):
        self.find_gravity()
        for body in self.bodies:
            body.update(self.dt, vel_change = self.check_collision(body))
#END of PhysEngine class
