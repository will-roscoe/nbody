
import nbody as nb
from nbody.tools.formatter import Formatter

# objects to query
objs = ('10','199','299','399','499','599','699','799','899')

# create a set of bodies at 24h later time
bphys = nb.Engine(dt=1800)
test_bodies = nb.horizons_batch(search_queries =objs , time = '2023-11-04')
bphys.attach_bodies(test_bodies)
nb.export_obj(bphys, 'data\\control_data')




# create a set of bodies and color for simulation later
bodies = nb.horizons_batch(search_queries =objs , time = '2023-11-03')
for i,color in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E',
                          '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    bodies[i].color = color
phys = nb.Engine(dt=1800)
phys.attach_bodies(bodies)
#choose only to simulate gravity between bodies
phys.do_bodygravity = True
phys.do_collisions = False
phys.do_fieldgravity = False
# export initial
nb.export_obj(phys, 'data\\initial_data')
#simulate 1 day
phys.simulate(48)
# export final
nb.export_obj(phys, 'data\\final_data')
diff_pos = [tb.pos-b.pos for (tb,b) in  zip(test_bodies, bodies)]
print([p-diff_pos[0] for p in diff_pos])
'''
RETURNED
(734208.3934936523, 1991346.5834197998, 95390.21323432401),
(30879.216217041016, -888387.130279541, -13982.13922296092),
(-388433.34998279624, -196322.83775445074, 10607.861980419606),
(127065.63578796387, 141509.94929504395, -151.14967669174075),
(-137375.16325378418, -43945.42427062988, -2756.059154611081),
(-19406.46745300293, 11261.982528686523, -3151.591775994748),
(-5973.261642456055, -433.52235412597656, -4189.279098611325),
(-3577.7086639404297, -2188.288528442383, 506.47766675427556)



'''







#save for visualisation
#phys.save_as('bodies', 'solarsys_bodies')
#phys.load_as('bodies', 'solarsys_bodies')


# test output periods
syst = nb.Engine(dt=1800)
bodies2 = nb.horizons_batch(search_queries =objs)
for i,col in enumerate(['#F28322','#BFBEBD', '#D9B391', '#63BAA6', '#F27A5E',
                          '#BFAE99', '#D9B779', '#95BBBF', '#789EBF']):
    bodies2[i].color = col 
syst.attach_bodies(bodies2)
syst.do_bodygravity = True
syst.do_collisions = False
syst.do_fieldgravity = False
syst.make_relative_to(bodies2[0])
syst.simulate(15000)

fmt = Formatter(items=['period'],engine=syst,plotskip=1,c_mass=syst.bodies[0].mass.c())
for b in syst.bodies:
    fmt.target = [b,100]
    print(str(fmt))
'''
output
Period: 0.2060 d (sun)
Period: 0.3234 a
Period: 0.6079 a
Period: 0.9765 a (earth)
Period: 1.9560 a
Period: 11.0644 a
Period: 30.4301 a
Period: 86.8747 a
Period: 163.4680 a
'''
dark = {'line':'w', 'face':(0,0,0,0), 'bkgd':'black', 'text':'w'}
sim = nb.MPLVisual(engine=syst, 
                name='Solar System', show_info=True,
                step_skip_frames=300, step_skip_points=150, max_period=1, do_picking=True, color_dict=dark,
                focus_body=syst.bodies[0])
sim.start()
