from mpl_toolkits.mplot3d.art3d import Text3D
from matplotlib.spines import Spine
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.axis3d import XAxis, YAxis, ZAxis
from .core import mplVisual
axes_3Ds = (XAxis, YAxis, ZAxis)
class Style:
    def __init__(self,gui:mplVisual, line_color, background_color, grid):
        gui.fig.subplots_adjust(0,0,1,1,0,0)
        if grid['show'] == False:
            self.ax.grid(False)
            ax_color = (0,0,0,0)
        else:
            ax_color = line_color
        
        for artist in (gui.fig, *gui.fig.get_children()):
            artist.set_facecolor(background_color)
        
        for artist in gui.ax.get_children():
            if isinstance(artist, Spine):
                artist.set_color(ax_color)
            if isinstance(artist, axes_3Ds):
                artist.set_pane_color((0.,0.,0.,0.))
            print(f'{artist}:{artist.get_children()}')
        
 
            #mpl.rcParams['text.color'] = self.pcolor
            

    def _axcolor(self, color):
        self.ax.tick_params(color=color,labelcolor=color)
        for ax in (self.ax.xaxis,self.ax.yaxis, self.ax.zaxis):
            ax.label.set_color(color)
        mpl.rcParams['axes.labelcolor'] = color
