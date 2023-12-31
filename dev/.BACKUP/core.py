import re
from datetime import datetime, timedelta
from decimal import Decimal
import math

from tqdm import tqdm, trange
import numpy as np

# Plotting and animation Packages
import matplotlib as mpl
#mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from cycler import cycler
from matplotlib.text import Text
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3D
# JPL Querying Packages
from astroquery.jplhorizons import Horizons
import astropy.units as u


from . import errors as e
from .base import (Variable, Vector, HistoricVector, _V, _O, typecheck,
Iterable, NoneType, NumType, VectorType, VarType, DecType)

try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)

def sphere(pos, radius, N=20):
    (c, r) = (pos, radius)
    u, v = np.mgrid[0:2*np.pi:N*1j, 0:np.pi:N*1j]
    x = r*np.cos(u)*np.sin(v) + c[0]
    y = r*np.sin(u)*np.sin(v) + c[1]
    z = r*np.cos(v) + c[2]
    return x,y,z
class Body:
    def __init__(self,
                mass,
                init_pos,
                init_vel=(0,0,0),
                radius = 0,
                bounce = 0.999,
                color= None,
                identity = None):
        if isinstance(identity, str):
            self.identity = identity
        elif identity is None:
            self.identity = f'body[m={mass}]'
        else:
            e.raise_type_error('identity', str, identity)
        
        init_pos, init_vel = _O(init_pos), _O(init_vel)
        typecheck(((init_pos,Iterable),(init_vel,Iterable),(mass,NumType),
                     (radius,NumType),(bounce,NumType),(color,(NoneType, str))))
        _id = self.identity
        self.pos = HistoricVector(li=init_pos, identity=f'{_id}_pos', units_v='m')
        self.vel = HistoricVector(li=init_vel, identity=f'{_id}_vel', units_v='ms^-1')
        self.acc = HistoricVector(0,0,0,identity=f'{_id}_acc',units_v='ms^-2')
        self.mass = Variable(mass,identity=f'{_id}_mass',units='kg')
        self.radius = Variable(radius,identity=f'{_id}_rad',units='m')
        self.bounce = bounce
        self.color = color
    def __str__(self):
        return f'Body("{self.identity}",\n\
    mass="{self.mass.c()} {self.mass.units}",\n\
    currentpos="{self.pos.c()} {self.pos.units}",\n\
    currentvel="{self.vel.c()} {self.vel.units}",\n\
    currentvel="{self.acc.c()} {self.acc.units}")'
    def __repr__(self):
        return f'Body("{self.identity}", m={self.mass.c()}, r={self.pos.c()},\
v={self.vel.c()}), a={self.acc.c()})'

    def update(self, dt = 1,
                vel_change=None,
                acc_change=None,
                vel_next=None):
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
            acc_change = _O(acc_change)
            if isinstance(vel_change, Iterable):
                if len(vel_change) == 3:
                    self.pos.next(self.pos + (Vector(vel) + vel_change)*dt)
                else:
                    raise RuntimeError
            else:
                raise RuntimeError
        else:
            self.pos.next(self.pos + Vector(vel)*dt)
        
    
    def _reinitialise(self, init_pos=None, init_vel=None):
        self.acc = HistoricVector(0,0,0,
                             identity=f'{self.identity}_acc',
                             units_v='ms^-2')
        if init_pos != None:
            if isinstance(init_pos, (*Iterable, *VectorType)):
                self.pos = HistoricVector(li=init_pos,
                                    identity=f'{self.identity}_pos',
                                    units_v='m')
            else:
                e.raise_type_error('init_pos', Iterable, init_pos)
        if init_vel != None:
            if isinstance(init_vel, (*Iterable, *VectorType)):
                self.vel = HistoricVector(li=init_vel,
                                    identity=f'{self.identity}_vel',
                                    units_v='ms^-1')
            else:
                e.raise_type_error('init_vel', (*Iterable, *VectorType), init_vel)

#END of Body Class

def horizons_query(searchquery,
                   observer = '0',
                   time = '2023-11-03',
                   num_type: type = float,
                   return_type = 'body'):
    typecheck(((searchquery, str),(observer,str), (time, str), (num_type, type), (return_type,str)))
    if all([return_type.lower() != opt for opt in ('body', 'dict', 'print')]) :
        e.raise_value_error('return_type', type, return_type)
    tqdm.write(f'Querying "{searchquery}" @JPL Horizons System')
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
    _raw0 = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors(get_raw_response=True)
    _tab = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors()
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _raw0)
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        tqdm.write(f'could not find name for object "{searchquery}", reverting to placeholder name.')
        name = f'JPL_{searchquery}'
    try:
        _raw =_raw0.split('Ephemeris')[0]
        m_exp_unit = _raw.split('Mass')[1].split('^')[1].split('=')[0].strip(') ,~ (')
        m_exp = "".join(c for c in m_exp_unit if c.isdigit() or c == '.')
        m_unit = "".join(c for c in m_exp_unit.lower() if c == 'k' or c == 'g')
        mass = "".join(c for c in  _raw.split('Mass')[1].split('^')[1].split('=')[1].strip(') ,~ (').lower() if 
                        c.isdigit() or c == '.' or c=='+' or c=='-')
    except IndexError:
        print(_raw0)
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

def horizons_batch(search_queries,
                   observer = '0',
                   time = '2023-11-03',
                   num_type= float,
                   return_type = 'body'):
    new_bodies = []
    for query in tqdm(search_queries, desc='Getting data from JPL Horizons', unit='queries'):
        if isinstance(query, str):
            new_bodies.append(horizons_query(query, observer, time, num_type, return_type))
        else:
            raise e.raise_type_error('item in search_queries', str, query)
    return new_bodies

#START of PhysEngine class
class PhysEngine:
    def __init__(self, dt = 1, checking_range=3):
        self.bodies, self.planes, self.fields = [], [], []
        typecheck(((dt, NumType), (checking_range,NumType)))
        self.dt = dt
        self._rangechk = checking_range
        self.do_collisions = True
        self.do_bodygravity = True
        self.do_fieldgravity = True
    def attach_bodies(self, new_bodies):
        if isinstance(new_bodies, Iterable):
            for i, new_body in enumerate(new_bodies):
                if isinstance(new_body, Body) or issubclass(type(new_body), Body):
                    self.bodies.append(new_body)
                else:
                    e.raise_type_error(f'new_body at index {i}', type(Body), new_body)
        else:
            e.raise_type_error('new_bodies', Iterable, new_bodies)
        tqdm.write(f'{len(self.bodies)} bodies currently attached to engine.')
    
    def _loadeng(self,eng):
        self = eng
    
    def save(self,dump='engine',file_name='nbody_data'):
        _saveobjs = {'bodies': {'bodies':self.bodies},
                     'engine':{'engine':self}}
        with open(f'{file_name}.npz','wb') as file:
            np.savez(file, **_saveobjs[dump])

    def load(self, objects='engine', file_name='nbody_data'):   
        _loadobjs = {'bodies': ("self.attach_bodies(objs['bodies'])",),
                     'engine': ("self._loadeng(objs['engine'])",)}
        with open(f'{file_name}.npz', 'rb') as file:
            objs = np.load(file, allow_pickle=True)
            for func in _loadobjs[objects]:
                try:
                    eval(func)
                except KeyError:
                    raise LookupError(f'cannot find {objects} value in "{file}"')

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
        _acc = _O(accel_vector)
        if len(_acc) == 3 and isinstance(_acc[0], NumType):
            self.fields+= [Vector(li=_acc)] 
        else:
            e.raise_type_error('accel_vector', (*Iterable, *VectorType), accel_vector)
    
    def create_plane(self, const_axis ='z', const_val  = 0):
        if const_axis in ('x', 'y', 'z') and isinstance(const_val, NumType):
            self.planes.append([const_axis, const_val])
            tqdm.write(f'constant plane {const_axis}={const_val} has been initialized.')
        else:
            e.raise_value_error('const_axis,const_val',(str,*NumType),(const_axis,const_val))
    
    def _check_collision(self, body, co_restitution=0):
        if self.do_collisions:
            returned_coll = False
            for bod in self.bodies:
                if body != bod: 
                    body_dist = (Vector(bod.pos.c()) - Vector(body.pos.c())).magnitude()
                                    
                    if body_dist < (bod.radius.c() + body.radius.c()):
                        returned_coll = True
                        n = Vector(bod.pos-body.pos)/body_dist
                        meff = 1/((1/bod.mass.c())+(1/body.mass.c()))
                        vimp = n*(body.vel - bod.vel)
                        imp = vimp*(1+co_restitution)*meff
                        dv = n*(-imp/body.mass.c())
                        return (dv+body.vel), False
            unitcomps = {'x':Vector((1,0,0)), 'y':Vector((0,1,0)), 'z':Vector((0,0,1))}
            for [pl_ax, pl_val] in self.planes:
                body_dists =  list(abs((Vector(body.pos.c()) +
                Vector(body.vel.c())*m*self.dt)[pl_ax] - pl_val) for m in range(self._rangechk))
                body_cur_dist = abs(body.pos[pl_ax] - pl_val)
                if any([(body_dists[i] < body.radius.c()) for i in range(self._rangechk)]) or\
                    body_dists[0] < body.radius.c():
                    if  body_dists[0]/body.radius.c() <= 1.01 and body_dists[-1]/body.radius.c() <= 1.01:
                        on_plane = True
                    else:
                        on_plane = False
                    pl_norm = unitcomps[pl_ax]
                    returned_coll = True
                    return (pl_norm*(-2*(body.vel*pl_norm)) + body.vel)*co_restitution, on_plane
            if not returned_coll:
                return Vector(body.vel.c()), False
        else:
            return Vector(body.vel.c()), False

    def _find_gravity(self, body):
        res = Vector((0,0,0))
        if self.do_bodygravity:
            for bod2 in self.bodies: 
                if bod2 != body:
                    dist_12 = body.pos - bod2.pos
                    unit_12, mag_12 = dist_12.unit(), dist_12.magnitude()
                    if isinstance(body.mass.c(), DecType) or isinstance(bod2.mass.c(), DecType):
                        b1m, b2m = Decimal(body.mass.c()), Decimal(bod2.mass.c())
                        force_on_bod1_t = (-1*Decimal(G) * b1m * b2m)
                        force_on_bod1 = unit_12*(force_on_bod1_t/Decimal(mag_12)**2)
                    else:
                        force_on_bod1_t = (-1*G * body.mass.c() * bod2.mass.c())
                        force_on_bod1 = unit_12*(force_on_bod1_t/(mag_12)**2)
                    res += force_on_bod1
        return res/body.mass.c()

    
    def evaluate(self):
        _temp = [[0,0,0] for _ in self.bodies]
        for i, body in enumerate(self.bodies):
            _temp[i][0:2] =self._check_collision(body, body.bounce)
        for i, body in enumerate(self.bodies):
            _temp[i][2] = self._find_gravity(body)
        for i, body in enumerate(self.bodies):
            col_vel, on_plane, acc_g = _temp[i]
            if not on_plane and self.do_fieldgravity:
                fieldvel = list(sum(f.c(i) for f in self.fields) for i in range(3))
            else:
                fieldvel = Vector((0,0,0)) 
            body.update(dt=self.dt, vel_next=(col_vel+fieldvel).c(), acc_change=acc_g)

# --- SUBCLASSES ---

# START of Simulation Class
class Simulation:
    def __init__(self,
                name = 'Nbody Simulation',
                engine = None,
                focus_body = None,
                focus_range = None,
                autoscale: bool = True,
                show_grid: bool = True, 
                show_shadows: bool = False,
                show_acceleration: bool = False,
                show_velocity: bool = False,
                vector_size = 1,
                labelling_type = 'legend',
                body_model = 'dots',
                guistyle = 'default', 
                do_picking = False, 
                info_body = None,
                show_info = False):
        typecheck(((engine,PhysEngine),(focus_range,(*NumType, NoneType)),
(autoscale,bool),(show_grid,bool),(show_shadows,bool),(show_acceleration,bool),(show_velocity,bool),
(vector_size,NumType),(labelling_type,str),(body_model,str),(guistyle,str), (do_picking, bool), (info_body, (NoneType, Body)),(show_info, bool), (focus_body, (NoneType, Body))))
        self.engine = engine
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
        self.focus_body = focus_body
        self.inf_b = info_body
        self.do_picking = do_picking
        self.show_info = show_info

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


    def _init(self, frames):
        for _ in trange(frames, desc='Evaluating motion for each frame', unit='frames'):
            self.engine.evaluate()
        #self.framelist = list(x for x in range(int(frames)))
        tqdm.write('Calculations finished, Starting interactive window...')
        print(self.engine.bodies, list(len(b.pos) for b in self.engine.bodies))
    
    def _animate(self, i):
        
        pr = 10
        co = {'clip_on':False}
        #maxim = len(self._engine.bodies[0].pos.X.record)
        ind = i#int(i*self._frameskip)
        step = 1#self._plotskip
        #while ind >= maxim:
           # ind =- 1 
        self.ax.clear()
        self.ax.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')
        
        if not self.autoscale:
            self.ax.set_box_aspect((1,1,1), zoom = self.zoom_slider.val)
            self.ax.set_autoscale_on(False)
            if self.focus_body is not None:
                (limx, limy, limz) = (float(m) for m in self.focus_body.pos[ind])
                if self.focus_range is None:
                    self.focus_range = float(max((max(self.focus_body.pos - bod.pos) for bod in self._engine.bodies)))
            else:
                limx, limy, limz = 0,0,0
                if self.focus_range is None:
                    self.focus_range = max(max(*bod.pos[ind]) for bod in self._engine.bodies)
            self.ax.set(xlim=((limx-self.focus_range),(limx+self.focus_range)),ylim=((limy-self.focus_range),
                        (limy+self.focus_range)),zlim=((limz-self.focus_range),(limz+self.focus_range)))
        else:
            self.ax.set_autoscale_on(True)
            self.ax.set_box_aspect(None, zoom = self.zoom_slider.val)
        
        xl, yl, zl, = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
        for plane in self.engine.planes:
            pl_const = np.array([[plane[1], plane[1]], [plane[1], plane[1]]])
            points = {'x':(pl_const,np.array([[yl[0],yl[1]],[yl[0],yl[1]]]),
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'y':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),pl_const,
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'z':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),
                           np.array([[yl[0],yl[0]],[yl[1],yl[1]]]),pl_const)}
            self.ax.plot_surface(*points[plane[0]], zorder=1,color=('xkcd:azure', 0.5), **co)
        
        for b in self.engine.bodies:
            _def = {'color':b.color, **co}
            _poshist = list(list(float(m) for m in _b.record[0:ind:step]) for _b in (b.pos.X, b.pos.Y, b.pos.Z))
            if self.show_velocity or self.show_acceleration:
                _pos = [float(m) for m in b.pos[ind]]
                if self.show_velocity:
                    _vel = [float(m) for m in b.vel[ind]]
                    self.ax.quiver(*_pos, *_vel, length=self.vector_size, color='red',
                                zorder=8, **co)
                if self.show_acceleration:
                    _acc = [float(m) for m in b.acc[ind]]
                    self.ax.quiver(*_pos, *_acc, length=self.vector_size, color='green',
                                zorder=8, **co)
            self.ax.plot(*_poshist, label=f'{b.identity}',
                         zorder=7, **_def)
            if self.body_model == 'dots':
                self.ax.scatter(*list(float(m) for m in b.pos[ind]),
                                marker='o', zorder=4, **_def)
            else:
                _sf = {'wireframe':self.ax.plot_wireframe, 'surface':self.ax.plot_surface}
                _sf[self.body_model](*sphere(list(float(m) for m in b.pos[ind]), b.radius.c()), zorder=2, **_def)

            if self.labelling_type == 'label':
                self.ax.text(*list(float(m) for m in b.pos[ind]), b.identity, zorder=10)
            if self.show_shadows:
                self.ax.plot(*_poshist[0:2],[(_poshist[2][ind]-self.focus_range)]*
                             len(_poshist[2]),
                    color='black', zorder=1.5, **co)
        if self.labelling_type == 'legend':
            self.leg = self.ax.legend(draggable=True)
            for text in self.leg.get_texts():
                text.set_picker(pr)
        def period(body):
            body = Body()
            m = 0
            for bod in self.engine.bodies:
                if bod.mass.c() > m:
                    m = bod.mass.c()
            a = max(body.pos[i].magnitude() for i in range(0,ind,self._plotskip))
            return 2*math.pi*math.sqrt((a**3)/(G*m))
        def ke(body):
            return (body.vel[ind])*(body.vel[ind])*(1/2)*body.mass.c() 
        self.ax.text2D(0.05, 0.05, f'ind:{ind}')
        if self.show_info and self.inf_b:
            self.ax.text2D(0.05, 0.95, f"{self.inf_b.identity}\npos:{self.inf_b.pos.c()}m\nvel:{self.inf_b.vel.c()}ms^-1\
\nacc:{self.inf_b.acc.c()}ms^-2\nmass:{self.inf_b.mass.c()} kg\nradius:{self.inf_b.radius.c()}\nKE:{ke(self.inf_b)} J\
\nperiod:{period(self.inf_b)} s", transform=self.inf_b)
    
    def _on_pick(self,event):
                if isinstance(event.artist,Line3D):
                    identity = event.artist.get_label()
                    self.focus_body = [b for b in self.engine.bodies if b.identity == identity]
                elif isinstance(event.artist, Text):
                    identity = event.artist.get_text()
                    self.inf_b = [b for b in self.engine.bodies if b.identity == identity]
                else:
                    return


    
    def start(self,eval_length=None,fps=None,frameskip=1,plotskip=1,cache=True):
        self._init(eval_length)
        self._plotskip, inv = plotskip, (1/fps)/1000
        tqdm.write('Starting Simulation Instance, Running Calculations:')
        anim = animation.FuncAnimation(self.fig,func=self._animate,init_func=self._init(eval_length), interval=inv,frames=eval_length,cache_frame_data=cache) 
        if self.do_picking:
            self.fig.canvas.mpl_connect('pick_event',self._on_pick)
        plt.show()
# --- SUBCLASSES ---
