import math
from ..tools import errors as e
from .base import typecheck, _O, _V, Iterable, NumType, NoneType, HistoricVector, Variable, Vector, VectorType

try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)

class Body:
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
            'radius':self.radius, 'color':self.color, 'bounce':self.bounce},'x':list(d.X for d in (self.pos, self.vel, self.acc)),
            'y':list(d.Y for d in (self.pos, self.vel, self.acc)),'z':list(d.Z for d in (self.pos, self.vel, self.acc)),
            '_data_':{'pos':self.pos['all'], 'vel': self.vel['all'], 'acc':self.acc['all'],'identity':self.identity, 'mass':self.mass, 
            'radius':self.radius, 'color':self.color, 'bounce':self.bounce}}[ind]

    def update(self,dt=1,vel_change=None,acc_change=None,vel_next=None):
        if vel_next:
            vel = _O(vel_next)
        else:
            vel = self.vel.c()
        
        if acc_change is not None:
            acc_change = _V(acc_change)
            self.vel.next((acc_change*dt + vel).c())
            self.acc.next(acc_change.c())

        else:
            self.vel.next((Vector((self.acc*dt))+vel).c())
            self.acc.next((0,0,0))
    
        if vel_change is not None:
            vel_change = _O(vel_change)
            
            if isinstance(vel_change,Iterable) and len(vel_change) == 3:
                self.pos.next(self.pos + (Vector(vel) + vel_change)*dt)
            else:
                raise RuntimeError
        else:
            self.pos.next(self.pos + Vector(vel)*dt)

    
    
    def _reinitialise(self,init_pos=None,init_vel=None):
        self.acc = HistoricVector(0,0,0,identity=f'{self.identity}_acc',units_v=self.acc.units)
        
        if init_pos != None:
            if isinstance(init_pos,(*Iterable,*VectorType)):
                self.pos = HistoricVector(li=init_pos,identity=f'{self.identity}_pos',units_v=self.pos.units)
            else:
                e.raise_type_error('init_pos',Iterable,init_pos)
        
        if init_vel != None:
            if isinstance(init_vel,(*Iterable,*VectorType)):
                self.vel = HistoricVector(li=init_vel,identity=f'{self.identity}_vel',units_v=self.vel.units)
            else:
                e.raise_type_error('init_vel',(*Iterable,*VectorType),init_vel)

    
    
    def get_(self, item, ind=-1, plotskip=0, c_mass=None, engine=None):
            def sma():
                if not plotskip >= ind:
                    a = max([(Vector(self.pos[i])-engine.barycenter(ind)).magnitude() for i in range(0,ind,plotskip)])
                    return a
                else:
                    return 'NaN'
            def per():
                a = sma()
                if a != 'NaN':
                    return 2*math.pi*math.sqrt((a**3)/(G*c_mass))
                else:
                    return a
            def ke():
                return Vector(self.vel[ind])*(self.vel[ind])*(1/2)*self.mass.c()
            
            _get_lookup = {**dict.fromkeys(['sma', 'semi_major_axis'],sma),
                       **dict.fromkeys(['per', 'period'],per),
                       **dict.fromkeys(['ke', 'kinetic_energy'], ke)}
            return _get_lookup[item]()
        