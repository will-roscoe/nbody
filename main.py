

import cProfile
import pstats
import nbody as nb
from line_profiler import profile
LINE_PROFILE=1

bodies = nb.horizons_batch(('10','199','299','399','499','599','699','799','899'))

for i,color in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E', '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    bodies[i].color = color
phys = nb.Engine(dt=1000)
phys.attach_bodies(bodies)
phys.make_relative_to(bodies[0])
phys.do_collisions = False
phys.do_fieldgravity = False
phys.simulate(20000)
#phys.save_as('bodies','solarsystem_bodies')


#phys.load_as('bodies','solarsystem_bodies')

sim = nb.mplVisual(engine=phys, 
                name='SS',
                focus_body=phys.bodies[0], show_info=True, autoscale=False, step_skip_frames=100, step_skip_points=100, max_period=3, cache=False, do_picking=True)


'''
import cProfile, pstats
profiler = cProfile.Profile()
profiler.enable()
main()
profiler.disable()
stats = pstats.Stats(profiler).sort_stats('ncalls').reverse_order()
stats.print_stats()

eng = nb.obj_from('./nbody/tests/testdata/test_enginerun.txt', 'engine')
print([bod['pos'] for bod in eng.bodies])
vis = nb.mplVisual(eng)
nb.export_obj(eng, 'data')
vis.start()
'''