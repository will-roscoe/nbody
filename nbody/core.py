
from decimal import Decimal

import math
from time import sleep
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
from matplotlib.text import Text



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
    currentvel="{self.acc.c()} {self.acc.units}")'
    
    
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
            'y':list(d.Y for d in (self.pos, self.vel, self.acc)),'z':list(d.Z for d in (self.pos, self.vel, self.acc))}[ind]

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

    
    
    def get_(self, item, ind, params):
            (plotskip, c_mass) = params
            def sma():
                try:
                    a = max(list(Vector(self.pos[ind]).magnitude() for i in range(0,ind,plotskip)))
                    return a
                except ValueError:
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
        




class Engine:
    def __init__(self,dt=1,checking_range=3):
        self.bodies, self.planes, self.fields= [], [], []
        (self.dt,self._rangechk) = typecheck(((dt,NumType),(checking_range,NumType)))
        self.do_collisions,self.do_bodygravity,self.do_fieldgravity = True,True,True
    
    def __len__(self):
        if len(self.bodies) is not None:
            return len(self.bodies[0].pos)
        else:
            return 0
    
    def attach_bodies(self, new_bodies):
        if isinstance(new_bodies,Iterable):
            for new_body in new_bodies:
                typecheck((new_body, Body))
                self.bodies.append(new_body)
            tqdm.write(f'{len(self.bodies)} bodies currently attached to engine.')
    
    
    def _loadeng(self,eng):
        self = eng
    
    
    def save_as(self,dump='engine',file_name='nbody_data'):
        _saveobjs = {'bodies':{'bodies':self.bodies},'engine':{'engine':self}}
        with open(f'{file_name}.npz','wb') as file:
            np.savez(file, **_saveobjs[dump])


    def load_as(self,objects='engine',file_name='nbody_data'):   
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

    
    def simulate(self,intervals):
        if isinstance(intervals, int) and len(self.bodies) > 0:
            for _ in trange(intervals, desc=f'Evaluating motion for each interval of {self.dt} seconds', unit='ints'):
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
            
            if intervals+1 == len(self.bodies[0].pos):
                tqdm.write(f'Finished Evaluating {intervals} intervals, ~{len(self.bodies[0].pos)} total intervals.')










PICKRADIUS = 10

_ax_loc = dict(zoom=(0.05,0.25,0.05,0.5), #ax1
               fps=(0.1,0.25,0.05,0.5), #ax2
               subadj=(0,0,1,1,0,0), #fig.subplots
               fsize=(16,9), b_as=(1,1,1)) #fig , fig.bbox

_c = dict(w='white', b='black', cl=(0.,0.,0.,0.),) 


_artists = dict(dot = dict(zorder=4, clip_on=False,picker=True, marker='o', pickradius=PICKRADIUS),
                plane=dict(zorder=1, clip_on=False, color=('xkcd:azure', 0.5)),
                trail=dict(zorder=7, clip_on=False, picker=True, pickradius=PICKRADIUS),
                surf=dict(zorder=2, clip_on=False, pickradius=PICKRADIUS),
                label=dict(zorder=10, clip_on=False),
                shadw=dict(zorder=1.5, clip_on=False, color='black'),
                vect=dict(zorder=8, clip_on=False))

_slider = dict(zoom=dict(label='Zoom', valmin=0.1, valmax=10, orientation='vertical'),
               fps=dict(label='Interval',valmax=2,orientation='vertical'))

_text = dict(info=dict(x=0.05, y=0.2, size='small'),
             ax_labels=dict(xlabel='x',ylabel='y',zlabel='z'),
             leg=dict(draggable=True, facecolor='black', fancybox=True))


class mplVisual:
    def __init__(self, engine, name='NBody Simulation (Matplotib)',
                focus_body=None,focus_range=None,
                autoscale=True,show_grid=True,show_shadows=False,
                show_acceleration=False,show_velocity=False,vector_size=1,
                labelling_type='legend',body_model='dots',guistyle='default',
                do_picking = False, show_info=False, info_params = {'output_raw':True,
                    'vector_pos':False, 'vector_vel':False, 'vector_acc':False}):
        
        (self.engine,self.focus_body,self.focus_range,self.autoscale,self.show_grid,
        self.show_shadows,self.show_acceleration,self.show_velocity,self.vector_size,
        self.labelling_type,self.body_model,self.guistyle,self.do_picking, self.show_info, self.info_params) = typecheck((
        (engine,Engine),(focus_body,(Body,NoneType)),(focus_range,(*NumType, NoneType)),
        (autoscale,bool),(show_grid,bool),(show_shadows,bool),(show_acceleration,bool),
        (show_velocity,bool),(vector_size,NumType),(labelling_type,str),(body_model,str),
        (guistyle,str),(do_picking,bool),(show_info, bool),(info_params, dict))) # new options here
        
        self.fig = plt.figure(name, figsize=_ax_loc['fsize'])
        self.ax = self.fig.add_subplot(projection='3d')
        self.ax1 = self.fig.add_axes(_ax_loc['zoom'])
        
        self.fig.subplots_adjust(*_ax_loc['subadj'])
        self.zoom_slider = Slider(self.ax1,valinit=1,**_slider['zoom'])
        (self._frameskip,self._plotskip, self.ax.computed_zorder, self.inf_b) = (
        1,1,False, self.focus_body)
        
        def _clearpanes():
            for ax in (self.ax.xaxis,self.ax.yaxis,self.ax.zaxis):
                ax.set_pane_color(_c['cl'])
        
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
                artist.set_facecolor(_c['b'])
            _clearpanes()
            _axcolor(_c['w'])
            mpl.rcParams['text.color'] = _c['w']
            self.tcolor = _c['w']
            self.zoom_slider.label.set_color(self.tcolor)
            
        else:
            self.tcolor = _c['b']
        
        if not self.show_grid:
            self.ax.grid(False)
            _axcolor(_c['cl'])
            _clearpanes()
        self._cmass = 0 # changed big mass thing here
        for bod in self.engine.bodies:
            if bod.mass.c() > self._cmass:
                self._cmass = bod.mass.c()
        if self.show_info:
            self.fm = Formatter(plotskip=1, c_mass=self._cmass, **self.info_params)
    
    def _draw_vectors(self,pos,other,c):
            self.ax.quiver(*pos,*other,length=self.vector_size,color=c,**_artists['vect'])

    def _animate(self,ind):
        sleep(self._int)
        self.ax.clear()
        
        # figure building

        self.ax.set(**_text['ax_labels'])
        if not self.autoscale:
            self.ax.set_box_aspect(_ax_loc['b_as'], zoom=self.zoom_slider.val)
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
        
        # plot planes

        for plane in self.engine.planes:
            xl, yl, zl, = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
            pl_const = np.array([[plane[1], plane[1]], [plane[1], plane[1]]])
            points = {'x':(pl_const,np.array([[yl[0],yl[1]],[yl[0],yl[1]]]),
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'y':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),pl_const,
                           np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                      'z':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),
                           np.array([[yl[0],yl[0]],[yl[1],yl[1]]]),pl_const)}
            self.ax.plot_surface(*points[plane[0]], **_artists['plane'])
            
        # plot bodies

        for b in self.engine.bodies:
            # cutting away any excessive looping orbits for performance
            cut = 0
            try:
                tau = (ind-(2*b.get_('period', ind, self._plotskip, self._cmass)/self.engine.dt))
                if tau > cut:
                    cut = math.ceil(tau)
            except TypeError:
                cut = 0
            # position and trail
            _poshist = list(list(float(m) for m in _b.record[cut:ind:self._plotskip]) for _b in (b.pos.X,b.pos.Y,b.pos.Z))
            _pos = [float(m) for m in b.pos[ind]]
            # draw vectors
            if self.show_velocity:
                self._draw_vectors(_pos,[float(m) for m in b.vel[ind]],'r')     
            if self.show_acceleration:
                self._draw_vectors(_pos,[float(m) for m in b.acc[ind]],'g')
            # draw trail
            self.ax.plot(*_poshist,label=f'{b.identity}',color=b.color, **_artists['trail'])
            # draw body
            if self.body_model == 'dots':
                self.ax.plot(*_pos,label=f'{b.identity}',color=b.color,**_artists['dot'])
            else:
                _sf = {'wireframe':self.ax.plot_wireframe, 'surface':self.ax.plot_surface}
                _sf[self.body_model](*sphere(_pos,b.radius.c()),label=f'{b.identity}', color=b.color **_artists['surf'])
            # draw label
            if self.labelling_type == 'label':
                self.ax.text(*_pos,b.identity,color=b.color, **_artists['label'])
            # draw shadows
            if self.show_shadows:
                self.ax.plot(*_poshist[0:2],[(_poshist[2][ind]-self.focus_range)]*
                             len(_poshist[2]), **_artists['shadw'])
        # draw legend
        if self.labelling_type == 'legend':
            handles, labels = self.ax.get_legend_handles_labels()
            handle_list, label_list = [], []
            for handle, label in zip(handles, labels):
                if label not in label_list:
                    handle_list.append(handle)
                    label_list.append(label)
            self.leg = self.ax.legend(handle_list, label_list, **_text['leg'])
            
            
            
            for text in self.leg.get_texts():
                text.set_picker(PICKRADIUS)
                text.set_color(self.tcolor)
        # draw info panel
        if self.show_info:
            self.fm.target = [self.inf_b, ind] # new all inside formatter object
            self.ax.text2D(s=str(self.fm), transform=self.ax.transAxes, **_text['info'])
    # object picking in interactive window
    def _on_pick(self,event):
            print(event.artist)
            if isinstance(event.artist,Text):
                identity = event.artist.get_text()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.inf_b = body
            elif isinstance(event.artist,Line3D):
                identity = event.artist.get_label()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.focus_body = body
            
  
    def start(self,interval=0.033,frameskip=1,plotskip=1,cache=False, speed_control=False):
        self._plotskip, self.int = plotskip, interval
        self.fm.par['ps'] = self._plotskip
        f = list(frameskip*x for x in range(int(len(self.engine)/frameskip)))
        tqdm.write('Starting Visual Environment')
        self._int = 0    
        anim = animation.FuncAnimation(self.fig, func=self._animate, interval=self.int*1000, frames=f, cache_frame_data=cache, ) 
        if speed_control:
            self.ax2 = self.fig.add_axes(_ax_loc['fps'])
            self.speed_slider = Slider(self.ax2,valinit=self.int, valmin=self.int,**_slider['fps'])
            self.speed_slider.label.set_color(self.tcolor)
            def _sp_ud(val):
                self._int = val-self.int
            self.speed_slider.on_changed(_sp_ud)
        if self.do_picking:
            self.fig.canvas.mpl_connect('pick_event',self._on_pick)
        plt.show()

from .text import Formatter