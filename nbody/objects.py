from body import Body

class Sol(Body):
    def __init__(self, init_pos: list | tuple = (((-1.232888890758686E+06)*(10**3)), ((-3.701017503718797E+05)*(10**3)), ((3.182189924211701E+04)*(10**3))),
                       init_vel: list | tuple = (((7.553670159837596E-03)*(10**3)), ((-1.314541679067443E-02)*(10**3)), ((-5.980735434994201E-05)*(10**3))))-> None:
        mass = 1989100000*(10**21)
        radius = 695508*(10**3)
        identity = 'Sol'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Jupiter(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3)), ((1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Saturn(Body):
    def __init__(self, init_pos: list | tuple = (((1.327777641282405E+09)*(10**3)), ((-6.014476320756447E+08)*(10**3)), ((-4.240744003900504E+07)*(10**3))),
                       init_vel: list | tuple = (((3.446779635869725E+00)*(10**3)), ((8.780998439535168E+00)*(10**3)), ((-2.903988774605155E-01)*(10**3))))-> None:
        mass = 5.6834*(10**26)
        radius = 60268*(10**3)
        identity = 'Saturn'
        super().__init__(mass, init_pos, init_vel, radius, identity)



'''
class Uranus(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Neptune(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Earth(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Venus(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Mars(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)

class Mercury(Body):
    def __init__(self, init_pos: list | tuple = (((5.674151411426915E+08)*(10**3)), ((4.794459035254424E+08)*(10**3)), ((-1.468375009239781E+07)*(10**3))),
                       init_vel: list | tuple = (((-8.577109696611455E+00)*(10**3)), ((1.059914774784837E+01)*(10**3), (1.479943675720525E-01)*(10**3))))-> None:
        mass = 1898187*(10**21)
        radius = 69911*(10**3)
        identity = 'Jupiter'
        super().__init__(mass, init_pos, init_vel, radius, identity)
'''
