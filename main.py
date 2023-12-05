

import nbody.core as nb
from nbody.errors import raise_type_error
import nbody.horizons as source
phys = nb.Engine(dt=1000)
'''
bodies = source.horizons_batch(('10','199','299','399','499','599','699','799','899'))

for i,color in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E', '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    bodies[i].color = color
phys.attach_bodies(bodies)
phys.make_relative_to(bodies[0])
phys.do_collisions = False
phys.do_fieldgravity = False
phys.simulate(20000)
phys.save_as('bodies','solarsystem_bodies')

'''
phys.load_as('bodies','solarsystem_bodies')

sim = nb.mplVisual(engine=phys, 
                   name='SS',
                   focus_body=phys.bodies[0], show_info=True, autoscale=False, frameskip=200, plotskip=200, max_p=1, max_pts=10, cache=False, do_picking=True) #do_picking=True, autoscale=False)

sim.start()

