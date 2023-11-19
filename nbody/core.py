import re
from datetime import datetime, timedelta
from multiprocessing import Pool, freeze_support
from decimal import Decimal
from matplotlib.text import Text
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
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3D


# JPL Querying Packages
from astroquery.jplhorizons import Horizons
import astropy.units as u


from . import errors as e
from .base import (NullVector, Variable, Vector, HistoricVector, _V, _O, typecheck,
Iterable, NoneType, NumType, VectorType, DecType)

try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)


def sphere(pos,radius,N=20):
    (c,r) = (pos,radius)
    u,v = np.mgrid[0:2*np.pi:N*1j, 0:np.pi:N*1j]

    x = r*np.cos(u)*np.sin(v)+c[0]
    y = r*np.sin(u)*np.sin(v)+c[1]
    z = r*np.cos(v)+c[2]
    
    return x,y,z








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
        
        self.pos = HistoricVector(li=init_pos,identity=f'{_id}_pos',units_v='m')
        self.vel = HistoricVector(li=init_vel,identity=f'{_id}_vel',units_v='ms^-1')
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
        self.acc = HistoricVector(0,0,0,identity=f'{self.identity}_acc',units_v='ms^-2')
        
        if init_pos != None:
            if isinstance(init_pos,(*Iterable,*VectorType)):
                self.pos = HistoricVector(li=init_pos,identity=f'{self.identity}_pos',units_v='m')
            else:
                e.raise_type_error('init_pos',Iterable,init_pos)
        
        if init_vel != None:
            if isinstance(init_vel,(*Iterable,*VectorType)):
                self.vel = HistoricVector(li=init_vel,identity=f'{self.identity}_vel',units_v='ms^-1')
            else:
                e.raise_type_error('init_vel',(*Iterable,*VectorType),init_vel)








def horizons_query(searchquery,observer='0',time='2023-11-03',num_type=float,return_type='body'):
    typecheck(((searchquery,str),(observer,str),(time,str),(num_type,type),(return_type,str)))
    if all([return_type.lower() != opt for opt in ('body','dict','print')]) :
        e.raise_value_error('return_type',type,return_type)
    
    tqdm.write(f'Querying "{searchquery}" @JPL Horizons System')
    _later_time = (datetime.strptime(time,'%Y-%m-%d')+timedelta(days=2)).strftime('%Y-%m-%d')
    
    _raw0 = Horizons(id=searchquery,location=observer,epochs={'start':time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors(get_raw_response=True)
    _tab = Horizons(id=searchquery,location=observer,epochs={'start':time,
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


def horizons_batch(search_queries,observer='0',time='2023-11-03',num_type=float,return_type='body'):
    new_bodies = []
    
    for query in tqdm(search_queries,desc='Getting data from JPL Horizons',unit='queries'):
        if isinstance(query,str):
            new_bodies.append(horizons_query(query,observer,time,num_type,return_type))
        else:
            raise e.raise_type_error('item in search_queries',str,query)
    
    return new_bodies

















class PhysEngine:
    def __init__(self,dt=1,checking_range=3):
        self.bodies = []
        (self.dt,self._rangechk) = typecheck(((dt,NumType),(checking_range,NumType)))
        self.do_collisions,self.do_bodygravity,self.do_fieldgravity = True,True,True
    
    
    def attach_bodies(self, new_bodies):
        if isinstance(new_bodies,Iterable):
            for new_body in new_bodies:
                typecheck((new_body, Body))
                self.bodies.append(new_body)
            tqdm.write(f'{len(self.bodies)} bodies currently attached to engine.')
    
    
    def _loadeng(self,eng):
        self = eng
    
    
    def save(self,dump='engine',file_name='nbody_data'):
        _saveobjs = {'bodies':{'bodies':self.bodies},'engine':{'engine':self}}
        with open(f'{file_name}.npz','wb') as file:
            np.savez(file, **_saveobjs[dump])


    def load(self,objects='engine',file_name='nbody_data'):   
        _loadobjs = {'bodies':("self.attach_bodies(objs['bodies'])",),
                     'engine':("self._loadeng(objs['engine'])",)}
        
        with open(f'{file_name}.npz','rb') as file:
            objs = np.load(file,allow_pickle=True)
            for func in _loadobjs[objects]:
                try:
                    eval(func)
                except KeyError:
                    raise LookupError(f'cannot find {objects} value in "{file}"')


    def make_relative_to(self,target_body):
        for body in self.bodies:
            if body != target_body:
                body._reinitialise((body.pos-target_body.pos), (body.vel-target_body.vel))
                body.identity = f'{body.identity}:Rel.{target_body.identity}'
        
        target_body._reinitialise((0,0,0),(0,0,0))
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

    def _check_collision(self,body,co_restitution=0):
        if self.do_collisions:
            returned_coll = False
            for bod in self.bodies:
                if body != bod: 
                    body_dist = (Vector(bod.pos.c()) - Vector(body.pos.c())).magnitude()             
                    if body_dist < (bod.radius.c() + body.radius.c()):
                        (returned_coll,n,meff) = (True,Vector(bod.pos-body.pos)/body_dist,1/((1/bod.mass.c())+(1/body.mass.c())))
                        dv = n*(-((n*(body.vel - bod.vel))*(1+co_restitution)*meff)/body.mass.c())
                        return (dv+body.vel), False
            
            unitcomps = {'x':Vector((1,0,0)),'y':Vector((0,1,0)),'z':Vector((0,0,1))}
            for pl in self.planes:
                body_dists = list(abs((Vector(body.pos.c()) +
                Vector(body.vel.c())*m*self.dt)[pl.c] - pl.v) for m in range(self._rangechk))
                if any([(body_dists[i] < body.radius.c()) for i in range(self._rangechk)]) or\
                body_dists[0] < body.radius.c():
                    on_plane = (True if body_dists[0]/body.radius.c() <= 1.01 and body_dists[-1]/body.radius.c() <= 1.01 else False)
                    returned_coll = True
                    return (unitcomps[pl.c]*(-2*(body.vel*unitcomps[pl.c])) + body.vel)*co_restitution, on_plane
            
            if not returned_coll:
                return Vector(body.vel.c()), False
        else:
            return Vector(body.vel.c()), False
    
    
    def _find_gravity(self,body):
        res = NullVector()
        if self.do_bodygravity:
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
                        force_on_bod1 = unit_12*(force_on_bod1_t/(mag_12)**2)
                    res += force_on_bod1
        
        return res/body.mass.c()

    
    def evaluate(self):
        _temp = [[0,0,0] for _ in self.bodies]
        
        for i,body in enumerate(self.bodies):
            _temp[i] = [*self._check_collision(body, body.bounce),self._find_gravity(body)]
        
        for i,body in enumerate(self.bodies):
            col_vel,on_plane,acc_g = _temp[i]
            if not on_plane and self.do_fieldgravity:
                fieldvel = list(sum(v.vec.c(i) for v in self.fields) for i in range(3))
            else:
                fieldvel = NullVector()
            body.update(dt=self.dt,vel_next=(col_vel+fieldvel).c(),acc_change=acc_g)




class Simulation:
    def __init__(self,name='Nbody Simulation',engine=None,
                focus_body=None,focus_range=None,
                autoscale=True,show_grid=True,show_shadows=False,
                show_acceleration=False,show_velocity=False,vector_size=1,
                labelling_type='legend',body_model='dots',guistyle='default',
                do_picking = False, show_info=False):
        
        (self.engine,self.focus_body,self.focus_range,self.autoscale,self.show_grid,
        self.show_shadows,self.show_acceleration,self.show_velocity,self.vector_size,
        self.labelling_type,self.body_model,self.guistyle,self.do_picking, self.show_info) = typecheck((
        (engine,PhysEngine),(focus_body,(Body,NoneType)),(focus_range,(*NumType, NoneType)),
        (autoscale,bool),(show_grid,bool),(show_shadows,bool),(show_acceleration,bool),
        (show_velocity,bool),(vector_size,NumType),(labelling_type,str),(body_model,str),
        (guistyle,str),(do_picking,bool), (show_info, bool)))
        
        self.fig = plt.figure(name, figsize=(16,9))
        self.ax = self.fig.add_subplot(projection='3d')
        self.ax1 = self.fig.add_axes((0.05,0.25,0.05,0.5))
        self.fig.subplots_adjust(left=0,bottom=0,right=1,top=1,wspace=0,hspace=0)
        (self.zoom_slider,self._frameskip,self._plotskip,self.ax.computed_zorder, self.inf_b) = (
        Slider(self.ax1,'Zoom',0.1,10,valinit=1,orientation='vertical'),1,1,False, self.focus_body)
        
        def _clearpanes():
            for ax in (self.ax.xaxis,self.ax.yaxis,self.ax.zaxis):
                ax.set_pane_color((0.,0.,0.,0.))
        
        def _axcolor(color):
            for spine in self.ax.get_children():
                if isinstance(spine,mpl.spines.Spine):
                    spine.set_color(color) 
            self.ax.tick_params(color=color,labelcolor=color)
            for ax in (self.ax.xaxis,self.ax.yaxis,self.ax.zaxis):
                ax.label.set_color(color)
            mpl.rcParams['axes.labelcolor'] = color
        
        if self.guistyle == 'dark':
            for artist in (self.fig,self.ax):
                artist.set_facecolor('black')
            _clearpanes()
            _axcolor('white')
            mpl.rcParams['text.color'] = 'white'
            self.tcolor = 'white'
        else:
            self.tcolor = 'black'
        
        if not self.show_grid:
            self.ax.grid(False)
            _axcolor((0.,0.,0.,0.))
            _clearpanes()
    
    
    def _init(self,eval_length,frameN,plotN):
        self.framelist = list(frameN*x for x in range(int(eval_length/frameN)))
        for _ in trange(eval_length, desc='Evaluating motion for each frame', unit='frames'):
            self.engine.evaluate()
        print(len(self.engine.bodies[0].pos))
        
        tqdm.write('Calculations finished, Starting interactive window...')
    
    
    def _draw_vectors(self,pos,other,c):
            self.ax.quiver(*pos,*other,length=self.vector_size,color=c,zorder=8,clip_on=False)
    

    
    def _animate(self,ind):
        co,pr = {'clip_on':False}, 10
        self.ax.clear()
        self.ax.set(xlabel='x',ylabel='y',zlabel='z')
        
        if not self.autoscale:
            self.ax.set_box_aspect((1,1,1), zoom=self.zoom_slider.val)
            self.ax.set_autoscale_on(False)
            
            if self.focus_body is not None:
                (limx,limy,limz) = (float(m) for m in self.focus_body.pos[ind])
                if self.focus_range is None:
                    self.focus_range = float(max((max(self.focus_body.pos - bod.pos) for bod in self.engine.bodies)))
            else:
                limx,limy,limz=0,0,0
                if self.focus_range is None:
                    self.focus_range = max(max(*bod.pos[ind]) for bod in self._engine.bodies)
            
            self.ax.set(xlim=((limx-self.focus_range),(limx+self.focus_range)),ylim=((limy-self.focus_range),
                        (limy+self.focus_range)),zlim=((limz-self.focus_range),(limz+self.focus_range)))
        else:
            self.ax.set_autoscale_on(True)
            self.ax.set_box_aspect(None,zoom=self.zoom_slider.val)
        
        #self.engine.const_objs.draw(self.ax)
        
        for b in self.engine.bodies:
            _def, _poshist = {'color':b.color, **co}, list(list(float(m) for m in _b.record[0:ind:self._plotskip]) for _b in (b.pos.X,b.pos.Y,b.pos.Z))
            _pos = [float(m) for m in b.pos[ind]]
            
            if self.show_velocity:
                self._draw_vectors(_pos,[float(m) for m in b.vel[ind]],'r')    
                
            if self.show_acceleration:
                self._draw_vectors(_pos,[float(m) for m in b.acc[ind]],'g')
            
            self.ax.plot(*_poshist,label=f'{b.identity}',zorder=7,**_def)
            
            if self.body_model == 'dots':
                self.ax.plot(*_pos,marker='o',zorder=4,pickradius= pr,**_def)
            else:
                _sf = {'wireframe':self.ax.plot_wireframe, 'surface':self.ax.plot_surface}
                _sf[self.body_model](*sphere(_pos,b.radius.c()),zorder=2,pickradius= pr,**_def)
            
            if self.labelling_type == 'label':
                self.ax.text(*_pos,b.identity,zorder=10)
            
            if self.show_shadows:
                self.ax.plot(*_poshist[0:2],[(_poshist[2][ind]-self.focus_range)]*
                             len(_poshist[2]),color='black',zorder=1.5,**co)
        
        def sm_a(body):
            try:
                a = max(list(Vector(body.pos[i]).magnitude() for i in range(0,ind,self._plotskip)))
                return a
            except ValueError:
                return 'NaN'
        
        def period(body):
            a = sm_a(body)
            if a != 'NaN':
                m = 0
                for bod in self.engine.bodies:
                    if bod.mass.c() > m:
                        m = bod.mass.c()
                return 2*math.pi*math.sqrt((a**3)/(G*m))
            else:
                return a
            
        
        
        def ke(body):
            return Vector(body.vel[ind])*(body.vel[ind])*(1/2)*body.mass.c()
        
        if self.labelling_type == 'legend':
            self.leg = self.ax.legend(draggable=True)
            for text in self.leg.get_texts():
                text.set_picker(pr)
        
        if self.show_info:
            self.ax.text2D(0.05, 0.2, f"{self.inf_b.identity}\npos:{self.inf_b.pos[ind]}m\nvel:{self.inf_b.vel[ind]}ms^-1\
\nacc:{self.inf_b.acc[ind]}ms^-2\nmass:{self.inf_b.mass.c()} kg\nbody radius:{self.inf_b.radius.c()}\nKE:{ke(self.inf_b)} J\
\nest period:{period(self.focus_body)} s", size='small', transform=self.ax.transAxes)

   
    def _on_pick(self,event):
            print(event.artist)
            if isinstance(event.artist,(Line3D, Text)):
                if isinstance(event.artist,Line3D):
                    identity = event.artist.get_label()
                else:
                    identity = event.artist.get_text()
                body = (b for b in self.engine.bodies if b.identity == identity)
                if not body:
                    return
                self.focus_body = body


    
    def start(self,eval_length=None,fps=None,frameskip=1,plotskip=1,cache=False):
        self._init(eval_length, frameskip, plotskip)
        self._plotskip, inv = plotskip, (1/fps)/1000
        tqdm.write('Starting Simulation Instance, Running Calculations:')
        anim = animation.FuncAnimation(self.fig,func = self._animate,interval=inv,frames=self.framelist,cache_frame_data=cache) 
        if self.do_picking:
            self.fig.canvas.mpl_connect('pick_event',self._on_pick)
        plt.show()

