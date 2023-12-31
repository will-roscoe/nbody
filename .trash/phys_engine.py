
try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)



from body import Body
from core import Vector, HistoricVector
import errors as e

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
    def __init__(self, dt: int | float = 1):
        self.bodies = []
        if isinstance(dt, (int, float)):
            self.dt = dt
        self.planes = []
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
        print(f'{len(self.bodies)} bodies attached')

    def create_plane(self, const_axis='z', const_val = 0):
        """
        Create a collision plane in the physics engine.

        Args:
            const_axis (str, optional): The constant axis for the collision plane ('x', 'y', or 'z').
            const_val (int or float, optional): The constant value for the collision plane.

        Raises:
            ValueError: If const_axis is not 'x', 'y', or 'z', or if const_val is not an int or float.
        """
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, (int, float)):
            self.planes.append([const_axis, const_val])
        else:
            e.raise_value_error('const_axis,const_val',((str),(int,float)),(const_axis,const_val))
    def check_collision(self,body):
        """
        Check for collisions between a body and other bodies or planes.

        Args:
            body (Body): The body to check for collisions.

        Returns:
            tuple: The resulting change in velocity for the body after collisions.

        Note:
            This method calculates and handles collisions with other bodies and collision planes.
        """
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
        """
        Calculate and apply gravitational forces between bodies.

        Note:
            This method calculates gravitational forces and updates the acceleration of each body.
        """
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
        """
        Evaluate the motion of attached bodies in the physics engine.

        Note:
            This method calculates the motion of attached bodies based on gravitational forces and collisions.
        """
        self.find_gravity()
        for body in self.bodies:
            body.update(self.dt, vel_change = self.check_collision(body))
#END of PhysEngine class
