
import nbody.core as nx
import nbody.models as nm
sim = nm.SolarSystemMB(show_info=True, do_picking=True)
sim.start(10000, frameskip=200, plotskip=200)