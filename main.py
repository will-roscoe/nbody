
import nbody as nb
'''
LINE_PROFILE=1

bodies = nb.horizons_batch(('10','199','299','399','499','599','699','799','899'))

for i,color in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E',
                          '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    bodies[i].color = color
phys = nb.Engine(dt=1000)
phys.attach_bodies(bodies)
phys.make_relative_to(bodies[0])
phys.do_collisions = False
phys.do_fieldgravity = False
phys.simulate(20000)
#phys.save_as('bodies','solarsystem_bodies')

'''
#phys.load_as('bodies','solarsystem_bodies')
eng= nb.Engine(dt=0.05)
bodies = (nb.Body(mass=5, init_pos=(-10.0,0,0), init_vel=(2.0,0,0), bounce=1, radius=2),
          nb.Body(mass=10, init_pos=(10.0,0,0), init_vel=(-2.0,0,0), bounce=1, radius=2),
            #Body(mass=1, init_pos=(15,1.,0), init_vel=(-0.5,0,0), bounce=1, radius=1)
            )
eng.do_bodygravity = False
eng.do_collisions = True
eng.do_fieldgravity = False
eng.attach_bodies(bodies)
eng.simulate(1000)
sim = nb.MPLVisual(engine=eng, 
                name='SS', show_info=True,
                step_skip_frames=1, step_skip_points=1, max_pts=1, cache=False, body_model='surface')

sim.start()
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