from core import HistoricVariable, HistoricVector, Vector
import errors as e


#START of Body Class
class Body:
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
