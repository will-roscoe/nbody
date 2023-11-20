import nbody.core as nb
import plotly.graph_objects as go
'''
bodies = nb.horizons_batch(('10','199','299','399','499','599','699','799','899'))
phys = nb.Engine(dt=1000)
phys.attach_bodies(bodies)
phys.make_relative_to(bodies[0])
phys.do_collisions = False
phys.simulate(100000)
phys.save_as('bodies')
'''
def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }
phys = nb.Engine(dt=1000)
phys.load_as('bodies')

(xs,ys,zs) = list(list(bod['pos'][x][::100] for bod in phys.bodies) for x in (0,1,2))

mm = [[0,0],[0,0],[0,0]]
for bod in phys.bodies[:5]:
    for i in (0,1,2):
        mm[i][0] = min(mm[i][0], *bod['pos'][i])
        mm[i][1] = max(mm[i][1], *bod['pos'][i])

b = phys.bodies
fig = go.Figure()
fb = b[0]
for i in range(5):        
    fig.add_trace(go.Scatter3d(x=[], y=[], z=[],mode="lines", line={'color':b[i].color}))
    #fig.add_trace(go.Scatter3d(x=[], y=[], z=[],mode="markers", marker=dict(color=b[i].color, size=10)))
    

fig.update(frames = [go.Frame(data= [go.Scatter3d(x=xs[i][:k+1],y=ys[i][:k+1],z=zs[i][:k+1])
                                     for i in range(5)],name=f'frame{k}') for k in range(1,len(xs[0]))])
#go.Scatter3d(x=[xs[i][k],],y=[ys[i][k],],z=[zs[i][k],]

fig.update_layout(scene=dict(
                    xaxis=dict(range=mm[0], autorange=False),yaxis=dict(range=mm[1], autorange=False),zaxis=dict(range=mm[2], autorange=False),
                    aspectratio=dict(x=1, y=1, z=1)),
      updatemenus=[dict(type="buttons",
                            buttons=[dict(label="Play",
                                            method="animate",
                                                args=[None, {'frame':{'duration':1}, 'transition':{'duration':1}}])])]
                                                )


fig.show()