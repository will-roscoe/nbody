from .core import PhysEngine, Simulation, Body, horizons_query

class SolarSystemMB(Simulation):
        def __init__(self, dt=1000,
                     engine = PhysEngine(),
                     show_grid= True,
                     show_shadows= False,
                     show_acceleration = False,
                     show_velocity= False,
                     vector_size = 1,
                     labelling_type = 'label',
                     body_model = 'dots',
                     guistyle = 'dark'):
            name = 'Major Bodies in Solar System'
            engine.dt = dt
            bodies = list(horizons_query(obj_id) for obj_id in (
                '10','199','299','399','499','599','699','799','899'))
            engine.attach_bodies(bodies)
            engine.make_relative_to(bodies[0])
            engine.do_collisions, engine.do_fieldgravity = False, False
            focus_body, focus_range, autoscale = bodies[0], None, False 
            super().__init__(name, engine, focus_body,
                             focus_range, autoscale, show_grid,
                             show_shadows, show_acceleration, show_velocity,
                             vector_size, labelling_type, body_model, guistyle)

        def start(self, eval_length=100000, fps=30, frameskip=5000, plotskip=5000):
            super().start(eval_length, fps, frameskip, plotskip)

class BouncingBalls(Simulation):
        def __init__(self, name = 'Balls in a Box',show_grid = True, show_shadows = False, show_acceleration = False, show_velocity = False, vector_size = 1, labelling_type = 'legend', guistyle = 'default'):
            balls = (
                    Body(3,(0,-1,1),(0,0.6,-0.01), 0.1, identity='4', color='yellow'),
                    Body(3,(0,1,0),(0.2,-0.5,-0.5), 0.1, identity='5', color='magenta'),
                    Body(3,(1,0,0),(-0.2,0.1,0.1), 0.1, identity='6', color = 'orange'),
                    Body(3,(0,-1,0),(-0.6,0.05,1), 0.1, identity='7', color='gold'))
            engine = PhysEngine()
            engine.do_collisions, engine.do_fieldgravity = True, True
            engine.dt = 0.01
            engine.attach_bodies(balls)
            for ax in ('x', 'y', 'z'):
                 for num in (-1.2, 1.2):
                    engine.create_plane(ax, num)  
            autoscale, body_model, focus_range, focus_body = False, 'surface', 1.2, None
            super().__init__(name, engine, focus_body, focus_range, autoscale, show_grid, show_shadows, show_acceleration, show_velocity, vector_size, labelling_type, body_model, guistyle)
            
        def start(self, eval_length=10000, fps=30, frameskip=10, plotskip=10):
            super().start(eval_length, fps, frameskip, plotskip)