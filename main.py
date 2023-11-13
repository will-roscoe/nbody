
import nbody.models as nb
bb = nb.BouncingBalls()
bb._engine.load('bodies')
print(bb._engine.bodies)
bb.start()