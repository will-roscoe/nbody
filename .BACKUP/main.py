
import nbody.core as nx

bodies = (nx.horizons_batch(('10', '399', '499')))
engine = nx.PhysEngine(dt=1000)
engine.do_collisions = False
engine.do_fieldgravity = False
engine.attach_bodies(bodies)
sim = nx.Simulation(engine=engine, focus_body=engine.bodies[0], guistyle='dark')

sim.start(100, 30, 1, 1)
