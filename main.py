
from doctest import BLANKLINE_MARKER
import random
import nbody.core as nb
def run_balls_in_box():
    balls = (
            nb.Body(3,(0,-1,1),(0,0.6,-0.01), 0.1, identity='4', color='yellow'),
            nb.Body(3,(0,1,0),(0.2,-0.5,-0.5), 0.1, identity='5', color='magenta'),
            nb.Body(3,(1,0,0),(-0.2,0.1,0.1), 0.1, identity='6', color = 'orange'),
            nb.Body(3,(0,-1,0),(-0.6,0.05,1), 0.1, identity='7', color='gold'))

    phys = nb.PhysEngine(dt=0.01)
    phys.attach_bodies(balls)
    phys.create_plane('z', -1.2)
    phys.create_plane('z', 1.2)
    phys.create_plane('x', -1.2)
    phys.create_plane('x', 1.2)
    phys.create_plane('y', -1.2)
    phys.create_plane('y', 1.2)

    sim = nb.Simulation('bouncing balls', phys, labelling_type='legend',
                        body_model='surface', autoscale=False, show_acceleration=True,
                        focus_range=1.2, vector_size=200)
    sim.start(frames=3000, fps=30, frameskip=2)

run_balls_in_box()