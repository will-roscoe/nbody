
import nbody.core as nb

phys = nb.Engine(dt=1000)
phys.load_as('bodies')
for i,color in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E', '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    phys.bodies[i].color = color
phys.save_as('bodies', 'solarsystem_bodies')
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
sim.start(fps=30, frameskip=1000, plotskip=200)