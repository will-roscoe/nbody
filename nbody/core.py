
from decimal import Decimal

import math
from time import sleep
from matplotlib.lines import Line2D
from tqdm import tqdm, trange
import numpy as np
import os

# Plotting and animation Packages
import matplotlib as mpl
#mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button, Slider
from cycler import cycler
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3D
from matplotlib.text import Text
from matplotlib.patches import Circle


from . import errors as e
from .base import (NullVector, Variable, Vector, HistoricVector, _V, _O, typecheck,
Iterable, NoneType, NumType, VectorType, DecType)


try:
    from scipy.constants import G
except ModuleNotFoundError:
    G = 6.6743*10**(-11)

PICKRADIUS = 10

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
            tqdm.write(f'«Engine» → {len(self.bodies)} bodies currently attached to engine.')
    
    
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
        tqdm.write(f"«Engine» → Bodies' positions and velocities have been made relative to {target_body.identity}.")
        target_body.identity = f'{target_body.identity}(Static)'
  

    
    def create_acceleration(self, accel_vector):
        _acc = _O(accel_vector)
        if len(_acc) == 3 and isinstance(_acc[0], NumType):
            self.fields+= [Vector(li=_acc)] 
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
                        force_on_bod1 = unit_12*(force_on_bod1_t/(mag_12)**2)
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
                        fieldvel = list(sum(v.vec.c(i) for v in self.fields) for i in range(3))
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
# kwargs
'''
kwargs

info_calc
show_acc
show_vel
vector_size
speed_control
start_ind


as dict
linecolor backgroundcolor facecolor textcolor

'''
        

class mplVisual:
    def __init__(self, 
                 engine,
                 name='NBody Simulation (Matplotib)',
                 show_info=False,
                 show_grid=True,
                 focus_body=None,
                 do_picking = True, 
                 **kwargs):
        self.args = {
            'color_dict': dict(line='black', face=(0,0,0,0), bkgd='white', text='black'),
            'speed_control': False,
            'vect_params': dict(vel=False, acc=False, size=1),
            'start_index':0,
            'info_calc':False,
            'focus_range':None,
            'anim_cache':False,
            'max_pts':None,
            'is_running':True, 
            'labelling_type':'legend',
            'body_model':'dots',
            'info_body':focus_body,
            'fmt_params':dict(),
            'file':None,
            'step_skip_frames':1,
            'step_skip_points':1,
            'max_period':2,
            'fps':30
            }
        self.files = dict()
        self.engine = engine
        self.show_info=show_info
        self.show_grid=show_grid
        self.focus_body=focus_body
        self.do_pick=do_picking
        self.args.update(kwargs)
        if self.show_grid == False:
            self.args['color_dict']['line'] = (0,0,0,0)
            self.args['color_dict']['face'] = (0,0,0,0)
        
        
        
        
        
        self.plt = dict(ptstep=self.args['step_skip_points'],
                        maxperiod=self.args['max_period'],
                        maxpts=self.args['max_pts'],
                        interval=1000/self.args['fps'],
                        frmstep=self.args['step_skip_frames'],
                        major_body=max((b.mass.c() for b in self.engine.bodies)))
        flist = list(self.plt['frmstep']*x for x in range(int(len(self.engine)/self.plt['frmstep'])))
        self.anim_args = dict(interval=(5 if self.args['speed_control'] is True else 1000/self.args['fps']), frames=flist, cache_frame_data=self.args['anim_cache'])
        
        self.trail_data = dict()
        with tqdm(total = len(self.engine.bodies)*len(flist), desc='«mplVisual» → Building Trails', unit='items') as tbar:
            for b in self.engine.bodies:
                body_data = dict()
                for f in flist: 
                    try:
                        tau = (f-(self.plt['frmstep']*b.get_('period', f, self.plt['ptstep'], self.plt['major_body'], engine=self.engine)/self.engine.dt))
                        if tau > 0:
                            lower = math.ceil(tau)
                        else: raise TypeError
                    except TypeError:
                        lower = self.args['start_index']
                    body_data_f = list(list(float(m) for m in _b.record[lower:f:self.plt['ptstep']]) for _b in (b.pos.X,b.pos.Y,b.pos.Z))
                    if self.plt['maxpts'] is not None:
                        while any(len(i) > self.plt['maxpts'] for i in body_data_f):
                            for i in body_data_f:
                                i.pop(0)
                    body_data[f] = body_data_f
                    tbar.update(1)
                self.trail_data[b] = body_data
        
        if show_info == True:
            self.fmt = Formatter(output_raw=False,items=['identity','mass','radius','energy','period','pos','vel','acc', 'time'], vector_pos=False, vector_vel=False, vector_acc=False, engine=self.engine,plotskip=self.plt['ptstep'], c_mass=self.plt['major_body'])
        
        if self.args['info_calc'] is True:
            self.info_data = dict()
            with tqdm(total = len(self.engine.bodies)*len(flist), desc='«mplVisual» → Precomputing All Descriptions', unit='items') as pbar:
                for b in self.engine.bodies: 
                    body_data = dict()
                    for f in flist:
                        self.fmt.target = [b,f]
                        body_data[f] = str(self.fmt)
                        pbar.update(1)
                    self.info_data[b] = body_data
        else:
            self.info_data = None    
        
        self.fig = plt.figure(name, figsize=(16,9))
        
        self.ax = self.fig.add_subplot(projection='3d', computed_zorder = False)
        
        self.zoom_ax = self.fig.add_axes((0.05,0.25,0.05,0.5))
        self.zoom_slider = Slider(self.zoom_ax,valinit=1,label='Zoom',valmin=0.1,valmax=10, orientation='vertical')
        
        self.plpa_ax = self.fig.add_axes((0.48,0.05,0.04, 0.05))
        self.playpause = Button(self.plpa_ax, label=' ▶ ▐▐ ', hovercolor='white', color=(0.5,0.5,0.5,0.5))
        self.playpause.label.set(fontstretch=0.6, fontsize='large', color='black')
        self.playpause.on_clicked(self._toggle)
        
        if self.args['speed_control'] == True:
            self.spd_ax = self.fig.add_axes((0.1,0.25,0.05,0.5))
            self.speed_slider = Slider(self.spd_ax,valinit=1000/self.plt['interval'], valmax=1000/5,label='Target\nFPS',valmin=0.1,orientation='vertical')
            self.speed_slider.label.set_color(self.args['color_dict']['text'])
            def _sp_ud(val):
                self.plt['interval'] = 1000/val
            self.speed_slider.on_changed(_sp_ud)
        
        self.fig.subplots_adjust(0,0,1,1,0,0)
        self.ax.tick_params(color=self.args['color_dict']['line'],labelcolor=self.args['color_dict']['text'])
        for artist in (self.fig,self.ax):
            artist.set_facecolor(self.args['color_dict']['bkgd'])
        for spine in self.ax.get_children():
            if isinstance(spine,mpl.spines.Spine):
                spine.set_color(self.args['color_dict']['line'])      
        mpl.rcParams['text.color'] = self.args['color_dict']['text']
        for ax in (self.ax.xaxis,self.ax.yaxis, self.ax.zaxis):
            ax.label.set_color((self.args['color_dict']['text'] if show_grid is True else (0,0,0,0)))
            ax.set_pane_color(self.args['color_dict']['face'])
        if show_grid == False:
            self.ax.grid(False)
        self.ax.set_autoscale_on(False)
        if self.args['save_to'] is not None:
            with open(f'{self.args['save_to']}.npz', 'wb') as file:
                np.savez(file, mplVisual=self)
                tqdm.write(f'«mplVisual» → Saved instance to {self.args['save_to']}.npz')
    
    def _draw_info(self, ind):
        if self.info_data is None:
            self.fmt.target = [self.args['info_body'], ind]
            inf_string = str(self.fmt)
        else:
            inf_string = self.info_data[self.args['info_body']][ind]
        self.ax.text2D(s=inf_string, transform=self.ax.transAxes, x=0.05, y=0.2, size='small', horizontalalignment='left',
                verticalalignment='bottom', color=self.args['color_dict']['text'])  
    def _draw_legend(self):
        handles, labels = self.ax.get_legend_handles_labels()
        handle_list, label_list = [], []
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        self.leg = self.ax.legend(handle_list, label_list, draggable=True, facecolor=self.args['color_dict']['bkgd'], fancybox=True, loc=1, fontsize='small')
        for text in self.leg.get_texts():
            text.set_picker(PICKRADIUS)
            text.set_color(self.args['color_dict']['text'])

    def _draw_vectors(self,pos,other,c):
        self.ax.quiver(*pos,*other,length=self.args['vect_params']['size'],color=c,zorder=8, clip_on=False)


    
    def _animate(self,ind):
        if self.args['speed_control'] == True:
            sleep((self.plt['interval']-5)/1000)
        self.ax.clear()
    
        self.ax.set(xlabel='$x$',
                    ylabel='$y$',
                    zlabel='$z$')
        
        
        self.ax.set_box_aspect((1,1,1),zoom=self.zoom_slider.val)
        
    
        if self.focus_body is not None:
            (limx,limy,limz) = (float(m) for m in self.focus_body.pos[ind])
            if self.args['focus_range'] is None:
                self.args['focus_range'] = float(max((max(abs(x) for x in (self.focus_body.pos - bod.pos)) for bod in self.engine.bodies)))
        else:
            limx,limy,limz=0,0,0
            if self.args['focus_range'] is None:
                self.args['focus_range'] = max(max(*bod.pos[ind]) for bod in self.engine.bodies)
        rng = self.args['focus_range']
        self.ax.set(xlim=((limx-rng),(limx+rng)),
                ylim=((limy-rng),(limy+rng)),
                zlim=((limz-rng),(limz+rng)))
        
        
        
        
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
            self.ax.plot_surface(*points[plane[0]], zorder=1, clip_on=False, color=('xkcd:azure', 0.5))
            
        # plot bodies

        for b in self.engine.bodies:
            _poshist = self.trail_data[b][ind]
            _pos = [float(m) for m in b.pos[ind]]
            # draw vectors
            if self.args['vect_params']['vel'] == True:
                self._draw_vectors(_pos,
                                   [float(m) for m in b.vel[ind]],
                                   'r')     
            if self.args['vect_params']['acc'] == True:
                self._draw_vectors(_pos,
                                   [float(m) for m in b.acc[ind]],
                                   'g')
            # draw trail
            self.ax.plot(*_poshist,
                         label=f'{b.identity}',
                         color=b.color,
                         zorder=7,
                         clip_on=False,
                         picker=True,
                         pickradius=PICKRADIUS)
            # draw body
            if self.args['body_model'] == 'dots':
                self.ax.plot(*_pos,
                             label=f'{b.identity}',
                             color=b.color,
                             zorder=4,
                             clip_on=False,
                             picker=True,
                             marker='o',
                             pickradius=PICKRADIUS)
            else:
                _sf = {'wireframe':self.ax.plot_wireframe, 'surface':self.ax.plot_surface}
                _sf[self.args['body_model']](*sphere(_pos,b.radius.c()),
                                   label=f'{b.identity}',
                                   color=b.color,
                                   zorder=2,
                                   clip_on=False,
                                   pickradius=PICKRADIUS)
            # draw label
        if self.args['labelling_type'] == 'label':
            for b in self.engine.bodies:
                self.ax.text(*[float(m) for m in b.pos[ind]],b.identity,
                            color=b.color,
                            zorder=10,
                            clip_on=False)
        else:
            self._draw_legend()
        # draw info panel
        if self.show_info == True:
            self._draw_info(ind)


# object picking in interactive window
    def _on_pick(self,event):
            tqdm.write(f'«mplVisual» → [click] {event.artist}')
            if isinstance(event.artist,Text):
                identity = event.artist.get_text()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.args['info_body'] = body
            elif isinstance(event.artist,(Line3D, Line2D)):
                identity = event.artist.get_label()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.focus_body = body
                self.focus_range = None
    
    def _toggle(self, event):
            if self.args['is_running'] == True:
                self.anim.event_source.stop()
                self.args['is_running'] = False
            else:
                self.anim.event_source.start()
                self.args['is_running'] = True        
    

    
    def start(self, **viewparams):
        tqdm.write('«mplVisual» → Starting Visual Environment')
        self.anim = animation.FuncAnimation(self.fig, func=self._animate, **self.anim_args) 
        if viewparams:
            self.ax.view_init(**viewparams)
        if self.do_pick:
            self.fig.canvas.mpl_connect('pick_event',self._on_pick)
        plt.show()
    
from .text import Formatter















