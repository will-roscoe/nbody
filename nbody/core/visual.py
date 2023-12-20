

#### Python Builtins
import math
# ⤤ efficient ceil function
from time import sleep
# ⤤ used when speed_control is True to change interval
#### 3rd Party Libs/Packages
from tqdm import tqdm
# ⤤ progress bars
import numpy as np
# ⤤ mgrid used to create surfaces for spheres
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button, Slider
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3D
from matplotlib.text import Text
# ⤤ visual library, importing multiple classes mainly for type checks on user clicks 
#### Local Imports
from ..tools.formatter import Formatter
# ⤤ used to parse output for info display when show_info is enabled


PICKRADIUS = 10

def sphere(pos,radius,n=20):
    (c,r) = (pos,radius)
    # get a cubic mesh of points
    u,v = np.mgrid[0:2*np.pi:n*1j, 0:np.pi:n*1j]
    # compute x,y,z values from spherical mesh
    x = r*np.cos(u)*np.sin(v)+c[0]
    y = r*np.sin(u)*np.sin(v)+c[1]
    z = r*np.cos(v)+c[2]
    return x,y,z
   

class MPLVisual:
    '''Object that outputs a visual, animated plot of an `core.Engine` instance, and can be customized to users\
        preferences.

>[!NOTE]
The `core.Engine` instance should already have computed a simulation before insterting into a `mplVisual` instance. 

### Parameters
|Parameter| Required |Type| Description|
|---|---|---| ---|
|`engine` | ✓ | `core.Engine` | engine containing the system that is to be visualized |
|`name` | ✕ | `str` | identifying name of the visual instance. |
|`show_info` | ✕ | `bool` | whether to show a description of live information of the `info_body`\
    (Default: `focus_body`).Default: `False`.|
|`show_grid` | ✕ | `bool` | whether to draw the grid and grid faces of the plot. Default:`True`.|
|`focus_body` | ✕ | `core.Body` | body to keep in the center of the plot. If left as `None`, the plot will not center\
      on a body. Default:`None`|
|`do_picking` | ✕ | `bool` | whether to enable object picking, allowing a user to choose a `focus_body` by clicking on \
    the object, or `info_body` by clicking on the label or legend item.|
### `kwargs`
|Parameter|Type|Default| Description|
|---|---|---|---|
|`step_skip_frames` | `int` |`1` | the step to iterate over of frames. reduces the amount of frames in the final\
      visual, but speeds up the animation. |
|`step_skip_points` | `int` | `1`| the step to iterate over of points. reduces the amount of points in each trail. \
    eg: if `step_skip_frames=10`, and `step_skip_points=5`, each frame would contain 2 more points.| 
|`fps` | `int` | `30`|framerate of the animation.|
|`max_pts` | `int` |`None` |maximum amount of points in each trail to plot. will remove oldest points if `max_pts` \
    is exceeded.|
|`max_period` | `int` |`2` |maximum amount of periods for an object to draw until points are cut.|
|||||
|`vect_params` | `dict` |`{vel:False, acc:False, size:1}` | parameters to show acceleration and velocity vectors and\
      a size scaling factor.|
|`speed_control` | `bool` |`False` |whether to enable a slider in the visual to control the speed.|
|`color_dict` | `dict` |`{line:'black', face:(0,0,0,0), bkgd:'white', text:'black'}` | color pallette used by the \
    visual.|
|`focus_range` | `base.Numtype` | `None`|range to plot the grid around, given `focus_body` is not `None`. if `None`,\
      the visual will autoscale, based on the furthest object.|
|`labelling_type` | `str` | `'legend'`|whether to use labels on objects or a legend. either `legend` or `labels`.|
|`body_model` | `str` |`'dots'` |how to draw the objects. can be `surface` or `wireframe` for spherical surface or \
    `dots` as markers.|
|||||
|`info_body` | `core.Body` |`focus_body` |initial object to show information for.|
|`info_calc` | `bool` |`False` |whether to compute all info strings for objects for all frames prior to animation.|
|`anim_cache` | `bool` |`False` |see `matplotlib.animation.FuncAnimation` parameter `cache`.|
|||||
|`is_running` | `bool` |`True` |value used to toggle animation using play/pause button in visual.|
|||||
|`fmt_params` | `dict` | see docstring|dictionary sent to the `Formatter` object to deal with the info output.|
|`file` | `str` |`None` |where to save the `mplVisual` object. if `None` then the object is not saved.|
|`start_index` | `int` |`0` |index of bodies' data to begin animation at.|

### Usage
#### `start(**viewparams)`
   - function to start instances output window and animation. 
   - `viewparams`: parameters to pass to `axes.view_init()`. initial viewing parameters as `elev`, `azim`, `roll`.
    '''
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
            'prebuild':False,
            'focus_range':None,
            'anim_cache':False,
            'max_pts':None,
            'is_running':True, 
            'labelling_type':'legend',
            'body_model':'dots',
            'info_body':focus_body,
            'fmt_params':dict(output_raw=False,
                              items=['identity','mass','radius','energy','period','pos','vel','acc', 'time'],
                              vector_pos=False, vector_vel=False, vector_acc=False),
            'step_skip_frames':1,
            'step_skip_points':1,
            'max_period':2,
            'fps':30,
            'file':None,
            }
        self.files = dict()
        self.engine = engine
        self.show_info=show_info
        self.show_grid=show_grid
        self.focus_body=focus_body
        self.do_pick=do_picking
        # update args dep. on which kwargs were defined.
        self.args.update(kwargs)
        # make gridlines and faces transparent if not show grid
        if self.show_grid == False:
            self.args['color_dict']['line'] = (0,0,0,0)
            self.args['color_dict']['face'] = (0,0,0,0)
        
    
        
        
        # plotting values
        self.plt = dict(ptstep=self.args['step_skip_points'],
                        maxperiod=self.args['max_period'],
                        maxpts=self.args['max_pts'],
                        interval=1000/self.args['fps'],
                        frmstep=self.args['step_skip_frames'],
                        major_body=max((b.mass.c() for b in self.engine.bodies)))
        # list of data intervals to create frames for
        flist = list(self.plt['frmstep']*x for x in range(int(len(self.engine)/self.plt['frmstep'])))
        # organise args for animation function.
        self.anim_args = dict(interval=(5 if self.args['speed_control'] is True else 1000/self.args['fps']),
                              frames=flist, cache_frame_data=self.args['anim_cache'])
        # build data for trails
        self.data = {}
        for b in self.engine.bodies:
            self.data[b] = list(list(float(m) for m in _b.record) for _b in (b.pos.X,b.pos.Y,b.pos.Z))
        if self.args['prebuild']:
            self.trail_data = dict()
            with tqdm(total = len(self.engine.bodies)*len(flist),
                    desc='«mplVisual» → Building Trails', unit='items') as tbar:
                for b in self.engine.bodies:
                    body_data = dict()
                    for f in flist: 
                        body_data[f] = self.gen_trail(f,b)
                        tbar.update(1)
                    self.trail_data[b] = body_data
        else:
            self.trail_data = None
        # init a formatter to manage the info readout
        if show_info == True:
            self.fmt = Formatter(engine=self.engine,plotskip=self.plt['ptstep'],
                                 c_mass=self.plt['major_body'], **self.args['fmt_params'])
        # if true precalculate all info readouts for each frame and body.
        if self.args['info_calc'] is True:
            self.info_data = dict()
            with tqdm(total = len(self.engine.bodies)*len(flist),
                      desc='«mplVisual» → Precomputing All Descriptions', unit='items') as pbar:
                for b in self.engine.bodies: 
                    body_data = dict()
                    for f in flist:
                        self.fmt.target = [b,f]
                        body_data[f] = str(self.fmt)
                        pbar.update(1)
                    self.info_data[b] = body_data
        else:
            self.info_data = None    
        # Figure instance
        self.fig = plt.figure(name, figsize=(16,9))
        # axes3D instance
        self.ax = self.fig.add_subplot(projection='3d', computed_zorder = False)
        # axes instance to contain zoom slider
        self.zoom_ax = self.fig.add_axes((0.05,0.25,0.05,0.5))
        self.zoom_slider = Slider(self.zoom_ax,valinit=1,label='Zoom',valmin=0.1,valmax=10, orientation='vertical')
        # axes instance to contain control button/s
        self.plpa_ax = self.fig.add_axes((0.48,0.05,0.04, 0.05))
        self.playpause = Button(self.plpa_ax, label=' ▶ ▐▐ ', hovercolor='white', color=(0.5,0.5,0.5,0.5))
        self.playpause.label.set(fontstretch=0.6, fontsize='large', color='black')
        self.playpause.on_clicked(self._toggle)
        # axes to contain speed control slider
        if self.args['speed_control'] == True:
            self.spd_ax = self.fig.add_axes((0.1,0.25,0.05,0.5))
            self.speed_slider = Slider(self.spd_ax,valinit=1000/self.plt['interval'],
                                       valmax=1000/5,label='Target\nFPS',valmin=0.1,orientation='vertical')
            self.speed_slider.label.set_color(self.args['color_dict']['text'])
            def _sp_ud(val):
                self.plt['interval'] = 1000/val
            self.speed_slider.on_changed(_sp_ud)
        # adjust fig to fill entire window
        self.fig.subplots_adjust(0,0,1,1,0,0)
        # style the fig and children
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
        if self.args['file'] is not None:
            with open(f'{self.args['file']}.npz', 'wb') as file:
                np.save(file, MPLVisual=self)
                tqdm.write(f'«mplVisual» → Saved instance to {self.args['save_to']}.npz')
    def _draw_info(self, ind):
        if self.info_data is None:
            # get info from formatter live
            self.fmt.target = [self.args['info_body'], ind]
            inf_string = str(self.fmt)
        else:
            # or read from precalculated.
            inf_string = self.info_data[self.args['info_body']][ind]
        self.ax.text2D(s=inf_string, transform=self.ax.transAxes,
                       x=0.05, y=0.2, size='small', horizontalalignment='left',
                verticalalignment='bottom', color=self.args['color_dict']['text'])  
    def _draw_legend(self):
        handles, labels = self.ax.get_legend_handles_labels()
        handle_list, label_list = [], []
        # only drawing single item in key for each body
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        self.leg = self.ax.legend(handle_list,
                                  label_list,
                                  draggable=True,
                                  facecolor=self.args['color_dict']['bkgd'],
                                  fancybox=True, loc=1, fontsize='small')
        for text in self.leg.get_texts():
            text.set_picker(PICKRADIUS)
            text.set_color(self.args['color_dict']['text'])

    def _draw_vectors(self,pos,other,c):
        self.ax.quiver(*pos,*other,length=self.args['vect_params']['size'],color=c,zorder=8, clip_on=False)
    
    def gen_trail(self, f, b):
        try:
            tau = (f-(self.plt['frmstep']*b.get_('period', f, self.plt['ptstep'],
                    self.plt['major_body'],engine=self.engine)/self.engine.dt))
            if tau > 0:
                lower = math.ceil(tau)
            else: 
                raise TypeError
        except TypeError:
            lower = self.args['start_index']
        body_data_f = list(list(m for m in 
                _b[lower:f:self.plt['ptstep']]) 
                for _b in self.data[b])
        if self.plt['maxpts'] is not None:
            while any(len(i) > self.plt['maxpts'] for i in body_data_f):
                for i in body_data_f:
                    i.pop(0)
        return body_data_f

    
    def _animate(self,ind):
        # if true, sleep to slow down fps manually
        if self.args['speed_control'] == True:
            sleep((self.plt['interval']-5)/1000)
        self.ax.clear()
    
        self.ax.set(xlabel='$x$',
                    ylabel='$y$',
                    zlabel='$z$')
        
        # set current zoom
        self.ax.set_box_aspect((1,1,1),zoom=self.zoom_slider.val)
        
    
        if self.focus_body is not None:
            # set plot range to max distances from focus body
            (limx,limy,limz) = (float(m) for m in self.focus_body.pos[ind])
            if self.args['focus_range'] is None:
                self.args['focus_range'] = float(max(
                    (max(abs(x) for x in (self.focus_body.pos - bod.pos)) for bod in self.engine.bodies)))
        else:
            limx,limy,limz=0,0,0
            if self.args['focus_range'] is None:
                # find body furthest from origin and centre plot on origine.
                self.args['focus_range'] = float(max(max(*bod.pos[ind]) for bod in self.engine.bodies))
        rng = self.args['focus_range']
        self.ax.set(xlim=((limx-rng),(limx+rng)),
                ylim=((limy-rng),(limy+rng)),
                zlim=((limz-rng),(limz+rng)))
        
        
        
        
        # plot planes
    
        for plane in self.engine.planes:
            xl, yl, zl, = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
            pl_const = np.array([[plane[1], plane[1]], [plane[1], plane[1]]])
            # build identity matrix-like dict
            points = {'x':(pl_const,np.array([[yl[0],yl[1]],[yl[0],yl[1]]]),
                            np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                        'y':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),pl_const,
                            np.array([[zl[0],zl[0]],[zl[1],zl[1]]])),
                        'z':(np.array([[xl[0],xl[1]],[xl[0],xl[1]]]),
                            np.array([[yl[0],yl[0]],[yl[1],yl[1]]]),pl_const)}
            self.ax.plot_surface(*points[plane[0]], zorder=1, clip_on=False, color=('xkcd:azure', 0.5))
            
        # plot bodies

        for i,b in enumerate(self.engine.bodies):
            _poshist = (self.gen_trail(ind,b) if self.args['prebuild'] is False else self.trail_data[b][ind])
            _pos = [float(m) for m in b.pos[ind]]
            _color = (b.color if all((b.color == n for n in (None,'None'))) else f'C{i}')
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
                         color=_color,
                         zorder=7,
                         clip_on=False,
                         picker=True,
                         pickradius=PICKRADIUS)
            # draw body
            if self.args['body_model'] == 'dots':
                self.ax.plot(*_pos,
                             label=f'{b.identity}',
                             color=_color,
                             zorder=4,
                             clip_on=False,
                             picker=True,
                             marker='o',
                             pickradius=PICKRADIUS)
            else:
                _sf = {'wireframe':self.ax.plot_wireframe, 'surface':self.ax.plot_surface}
                _sf[self.args['body_model']](*sphere(_pos,b.radius.c()),
                                   label=f'{b.identity}',
                                   color=_color,
                                   zorder=2,
                                   clip_on=False,
                                   pickradius=PICKRADIUS)
            # draw label
        if self.args['labelling_type'] == 'label':
            for b in self.engine.bodies:
                self.ax.text(*[float(m) for m in b.pos[ind]],b.identity,
                            color=_color,
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
            # if text, update info readout
            if isinstance(event.artist,Text):
                identity = event.artist.get_text()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.args['info_body'] = body
            # if line/body update focus point
            elif isinstance(event.artist,(Line3D, Line2D)):
                identity = event.artist.get_label()
                body = list(b for b in self.engine.bodies if b.identity == identity)[0]
                if not body:
                    return
                self.focus_body = body
                self.focus_range = None
    
    def _toggle(self, event):
            # start/stop animation with control button
            if self.args['is_running'] == True:
                self.anim.event_source.stop()
                self.args['is_running'] = False
            else:
                self.anim.event_source.start()
                self.args['is_running'] = True        
    

    def start(self, **viewparams):
        ''' 
    - function to start instances output window and animation. 
   - `viewparams`: parameters to pass to `axes.view_init()`. initial viewing parameters as `elev`, `azim`, `roll`.'''
        tqdm.write('«mplVisual» → Starting Visual Environment')
        self.anim = animation.FuncAnimation(self.fig, func=self._animate, **self.anim_args) 
        # change view point if viewparams specified
        if viewparams:
            self.ax.view_init(**viewparams)
        # connect picking callback
        if self.do_pick:
            self.fig.canvas.mpl_connect('pick_event',self._on_pick)
        plt.show()
















