




import nbody.core as nb
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
sim = nb.mplVisual(phys, 'SS', phys.bodies[0],None, False,
                     show_grid= True,
                     show_shadows= False,
                     show_acceleration = False,
                     show_velocity= False,
                     vector_size = 1,
                     labelling_type = 'legend',
                     body_model = 'dots',
                     guistyle = 'dark',
                     do_picking = True,
                     show_info = True)
sim.start(frameskip=1000, plotskip=200, speed_control=True, cache=True)
