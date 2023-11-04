from core import Vector
from simulation import Simulation
from phys_engine import PhysEngine
from body import Body
from objects import *
import matplotlib.patches as mpl
bodies= (Saturn(), Sol(), Jupiter())
phys = PhysEngine(dt=1000000)
phys.attach_bodies((bodies))
sim = Simulation('solarsystemtest', engine=phys)
sim.start(duration=20, fps=90)
