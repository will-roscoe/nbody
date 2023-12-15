import pytest

from .. import Engine, Body
import mpmath as mp

class TestBody:
    def test_init_body(self):
        _body = Body(mass=1, init_pos=(0,0,0), init_vel=(1,0,0), radius=1)
        assert True
    def test_vel(self):
        body = Body(mass=1, init_pos=(0,0,0), init_vel=(1,0,0), radius=1)
        body.update(1)
        c_pos = body.pos.c()
        body.update(1, (1,0,0))
        nc_vel = body.vel.c()
        nc_pos = body.pos.c()
        assert c_pos == (1,0,0) and nc_pos == (3,0,0) and nc_vel == (1,0,0)

class TestEngine:
    def test_init_eng(self):
        _eng = Engine()
        assert True
    def test_eng_grav_f(self):
        eng = Engine()
        eng.attach_bodies((Body(mass=1, init_pos=(0,0,0), init_vel=(1,0,0), radius=1),Body(1000, (100,0,0), radius=10)))
        eng.do_bodygravity = True
        eng.do_collisions, eng.do_fieldgravity = False, False
        fg = eng._find_gravity(eng.bodies[0])
        assert mp.almosteq(fg.c(0), 6.674e-12)
        # fg almost equal to 6.674e-12 0, 0
    def test_eng_grav_f2(self):
        eng2=Engine()
        eng2.attach_bodies((Body(mass=1, init_pos=(0,0,0), init_vel=(1,0,0), radius=1),
                            Body(1000, (100,0,0), radius=10),
                            Body(mass=1000, init_pos=(-50, 0, 100), radius=10)))
        eng2.do_bodygravity = True
        eng2.do_collisions, eng2.do_fieldgravity = False, False
        fg2 = eng2._find_gravity(eng2.bodies[0])
        print(fg2)
        assert mp.almosteq(fg2.c(0), 4.2864e-12) and mp.almosteq(fg2.c(2), 4.7757e-12)
    def test_body_collisions(self):
        eng3 = Engine()
        bodies = (Body(mass=5, init_pos=(-10.0,0,0), init_vel=(2.0,0,0), bounce=1, radius=2),
                  Body(mass=10, init_pos=(10.0,0,0), init_vel=(-2.0,0,0), bounce=1, radius=2),
                  Body(mass=1, init_pos=(15,1.,0), init_vel=(-0.5,0,0), bounce=1, radius=1),
                 )
        eng3.do_bodygravity = False
        eng3.do_collisions = True
        eng3.do_fieldgravity = True
        eng3.attach_bodies(bodies)
        eng3.simulate(1000)
        _ke = list(sum(b.get_('ke',ind=i) for b in eng3.bodies) for i in (0,-1))
        assert mp.almosteq(_ke[0], _ke[1])