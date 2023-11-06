import errors as e
from core import NoneType
from body import Body
from phys_engine import PhysEngine

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from cycler import cycler
class Simulation:
    """
    A class for visualizing and controlling a physics simulation.

    This class allows you to create a 3D visualization of a physics simulation using
    Matplotlib and provides controls for adjusting the simulation parameters and viewing
    options. It includes features such as zoom, grid visibility, shadows, acceleration and
    velocity vectors, body models, and labeling options.

    Args:
        name (str): The name of the simulation.
        engine (PhysEngine|NoneType): The physics engine to simulate the physics of the system.
        focus_body (Body|NoneType): The body to focus on during the simulation.
        focus_range (int, float): The range for focusing on the body.
        autoscale (bool): Enable or disable autoscaling.
        show_grid (bool): Show or hide the grid in the visualization.
        show_shadows (bool): Show or hide shadows.
        show_acceleration (bool): Show or hide acceleration vectors.
        show_velocity (bool): Show or hide velocity vectors.
        vector_size (int, float): The size of velocity and acceleration vectors.
        labelling_type (str): The type of labeling, either 'legend' or 'label'.
        body_model (str): The model for representing bodies, such as 'dots', 'wireframe', or 'surface'.
        guistyle (str): The GUI style, either 'default' or 'dark'.

    Methods:
        start(frames, interval, duration, fps): Start the simulation animation.
    
    """
    def __init__(self,
                name: str = 'Nbody Simulation',
                engine: PhysEngine|NoneType = None,
                focus_body: Body|NoneType = None,
                focus_range: int|float = 0.5,
                autoscale: bool = True,
                show_grid: bool = True, 
                show_shadows: bool = False,
                show_acceleration: bool = False,
                show_velocity: bool = False,
                vector_size: int|float = 1,
                labelling_type: str = 'legend',
                body_model: str = 'dots',
                guistyle: str = 'default'):
        _argspec = {engine:(PhysEngine,NoneType),focus_range:(int,float),
autoscale:bool,show_grid:bool,show_shadows:bool,show_acceleration:bool,show_velocity:bool,vector_size:(int,float),
labelling_type:str,body_model:str,guistyle:str}
        
            
        for i,(arg,typ) in enumerate(_argspec.items()):
            if not isinstance(arg, typ):
                e.raise_type_error(f'Simulation arg({i})', typ, arg)
                break
        else:
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
        if issubclass(type(focus_body), Body) or isinstance(focus_body, (NoneType, Body)):
            self.focus_body = focus_body
        else:
            e.raise_type_error('focus_body', (Body, 'subclassBody', NoneType), focus_body)

        
        self.fig = plt.figure(name, figsize=(16,9))
        self.ax = self.fig.add_subplot(projection = '3d')
        self.ax1 = self.fig.add_axes((0.05,0.25,0.05,0.5))
        self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.zoom_slider = Slider(self.ax1, 'Zoom', 0.1, 10, valinit=1, orientation='vertical')
        
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


    def _init(self, i):
        for _ in range(i):
            self.engine.evaluate()

    def _animate(self, i):
        self.ax.clear()
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')
        self.ax.set_box_aspect((1,1,1), zoom = self.zoom_slider.val)
        if not self.autoscale:
            self.ax.set_aspect('equal', adjustable='box')
            self.ax.set_autoscale_on(False)
            if self.focus_range is not None and self.focus_body is not None:
                self.ax.set_xlim(xmin=(self.focus_body.pos.X[i]-self.focus_range),
                                 xmax=(self.focus_body.pos.X[i]+self.focus_range))
                self.ax.set_ylim(ymin=(self.focus_body.pos.Y[i]-self.focus_range),
                                 ymax=(self.focus_body.pos.Y[i]+self.focus_range))
                self.ax.set_zlim(zmin=(self.focus_body.pos.Z[i]-self.focus_range),
                                 zmax=(self.focus_body.pos.Z[i]+self.focus_range))
        else:
            self.ax.set_autoscale_on(True)
            self.ax.set_autoscalez_on(True)
        if self.show_velocity or self.show_acceleration:
            _pos = [[b.pos.X[i] for b in self.engine.bodies],
            [b.pos.Y[i] for b in self.engine.bodies],[b.pos.Z[i] for b in self.engine.bodies]]
            if self.show_velocity:
                _vel = [[b.vel.X[i] for b in self.engine.bodies],[b.vel.Y[i] for b in self.engine.bodies],
                [b.vel.Z[i] for b in self.engine.bodies]]
                self.ax.quiver(*_pos, *_vel, length=self.vector_size, color='red')
            if self.show_acceleration:
                _acc = [[b.acc.X[i] for b in self.engine.bodies],
                [b.acc.Y[i] for b in self.engine.bodies],[b.acc.Z[i] for b in self.engine.bodies]]
                self.ax.quiver(*_pos, *_acc, length=self.vector_size, color='green')
        for b in self.engine.bodies:
            self.ax.plot(b.pos.X.hist[0:i], b.pos.Y.hist[0:i], b.pos.Z.hist[0:i],
                        label=f'{b.identity}', zorder=2.9)
            if self.body_model == 'dots':
                self.ax.scatter(b.pos.X[i], b.pos.Y[i], b.pos.Z[i],
                                marker='o', zorder=4)
            if self.body_model == 'wireframe':
                pass
            if self.body_model == 'surface':
                pass
            if self.labelling_type == 'label':
                self.ax.text(*b.pos[i], b.identity)
            if self.show_shadows:
                self.ax.plot(b.pos.X.hist[0:i], b.pos.Y.hist[0:i],
                    [(self.focus_body.pos.Z[i]-self.focus_range)]*len(b.pos.Z.hist[0:i]),
                    color='black', zorder=2.9)
        if self.labelling_type == 'legend':
            self.ax.legend()


    def start(self, frames=None, interval=None, duration=None, fps=None):
        if frames and interval:
            f,i = frames,interval
        elif duration and fps:
            f,i = duration*fps, 1/fps 
        elif frames and fps:
            f,i = frames, 1/fps
        elif duration and interval:
            f,i = duration/interval, interval
        anim = animation.FuncAnimation(self.fig, func = self._animate, init_func = self._init(f), interval=i, frames=f) 
        plt.show()
        
