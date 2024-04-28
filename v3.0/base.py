
# numpy based vector objects, fixed dtype float64

# record vectors need to have variable length, but numpy won't allow this with np.array .
# we could use a list of numpy arrays but it may be slow to access elements.
# we need to be able to acess elements quickly within a numpy system with variable length, so 

import re
from typing import List, Union
import numpy as np
import random
from vispy import scene
from vispy.scene import visuals
class Particle:
    def __init__(self, pos: np.ndarray, vel: np.ndarray, mass: float, radius: float, charge: float, color: Union[str, np.ndarray] = 'white'):
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.mass = mass
        self.radius = radius
        self.charge = charge
        self.color = color
    def __repr__(self) -> str:
        return f"\nParticle(p={self.pos},v={self.vel},m={self.mass},r={self.radius},c={self.charge})"

class Engine:
    G = 6.67430e-11  # gravitational constant
    k = 8.9875517923e9  # Coulomb's constant
    mu0 = 4*np.pi*1e-7  # permeability of free space

    def __init__(self, particles: List[Particle], dt: float = 1e-3):
        self.particles = particles
        self.dt = dt
    def update(self):
        for p in self.particles:
            p.pos += p.vel * self.dt
            p.vel += self.force(p, self.dt) * self.dt / p.mass

    def resultantforce_nbody(self, p: Particle) -> np.ndarray: # slow method to calcluate resultant force on a particle
        force = np.zeros(3)
        for other in self.particles:
            if other is not p:
                r = other.pos - p.pos
                r_mag = np.linalg.norm(r)
                # gravitational force
                force += self.G * p.mass * other.mass / r_mag**3 * r
                # electric force
                force += self.k * p.charge * other.charge / r_mag**3 * r
                # magnetic force due to moving charges
                dB = self.mu0/(4*np.pi) * (other.charge * other.vel * self.dt * np.cross(r, other.vel)) / r_mag**3
                force += p.charge * np.cross(p.vel, dB)
        return force
class Node:
    def __init__(self, pos: np.ndarray, size: float, mass: float=0., depth=0) -> None:
        self.pos = pos
        self.size = size
        self.mass = mass
        self.children = []  # child nodes
        self.particles = []  # particles in this node
        self.is_end = False
        self.depth = depth
        self.cmass = pos
        #debug print(f"d{depth} Node: pos={self.pos}, size={self.size}, mass={self.mass}")
    def insert(self, particle: Particle) -> None:
        self.particles.append(particle)
        self.mass += particle.mass
        #debug print(f"inserted particle with mass {particle.mass} into node with mass {self.mass}")

    def split(self) -> None:
        if not self.is_end:
            #debug print(f"splitting a d{self.depth} Node: pos={self.pos}, size={self.size}, mass={self.mass}")
            
            #split into 8 octants
            for i in [1,-1]:  # bisecting in each dimension
                for j in [1,-1]:
                    for k in [1,-1]:
                        pos = self.pos + np.array([i, j, k]) * self.size / 2
                        self.children.append(Node(pos, self.size/2, depth=self.depth+1))  # add as child node
            # redistribute particles to children
            for particle in self.particles: 
                for child in self.children:  
                    if child.contains(particle):
                        child.insert(particle)
                        break
                else:
                    raise ValueError(f"particle not in any child node: \n{particle}")
                    
            self.children = [child for child in self.children if child.mass > 0]
            self.particles = []  # clear particles list after redistribution
            for child in self.children:
                if len(child.particles) > 1:
                    child.split()
                    self.fast_cmass()
                else:
                    child.is_end = True
    def contains(self, particle: Particle) -> bool:
        return np.all(np.abs(particle.pos - self.pos) <= self.size)   
    
    def fast_cmass(self) -> np.ndarray:
        self.cmass = np.sum([p.mass*p.pos for p in [*self.particles,*self.children]], axis=0) / self.mass
        return self.cmass
    def __repr__(self) -> str:
        return (f"""
d{self.depth}Node: pos={self.pos} size={self.size} mass={self.mass} 
children({len(self.children)}): {[child for child in self.children]}
particles({len(self.particles)}): {[particle for particle in self.particles]}""")

    def treeview(self) -> str:
        bars = "   ⎸"*self.depth
        info = f"⟶  d{self.depth} Node[{len(self.children)}ch,{len(self.particles)}p], cmass {self.cmass}, Mass:{self.mass}"
        if not self.is_end:
            info = "\u001b[37m\u001b[1m" + info + "\u001b[37m\u001b[0m"
        print(bars + info)
        for c in self.children:
            c.treeview()

    def dist(self, other):
        return np.linalg.norm(self.cmass - other.cmass)
class BHTree:
    def __init__(self, particles, origin:np.ndarray=np.zeros(3), size:float=100.) -> None:
        [np.mean(comp) for comp in zip(*[p.pos for p in particles])]
        #initpos = np.array([np.mean(comp) for comp in zip(*[p.pos for p in particles])])
        self.root = Node(origin, #initpos, 
                         size #np.max([np.max([np.abs(particle.pos - initpos)]) for particle in particles])
                         )
        for particle in particles:
            self.root.insert(particle)
        self.root.split()

        
        
    # barnes hut algorithm

    # 1, build octree by
        # defining a boundary volume
        # split volume into 8 octants
        # if octant contains >1 particles, split again, otherwise leave octant as leaf
        # continue until all octants contain 0 or 1 particles  


'''
Calculating the force acting on a body
To calculate the net force on a particular body, the nodes of the tree are traversed, starting from the root.
If the center of mass of an internal node is sufficiently far from the body, the bodies contained in that part
of the tree are treated as a single particle whose position and mass is respectively the center of mass and total
mass of the internal node. If the internal node is sufficiently close to the body, the process is repeated for
each of its children.

Whether a node is or isn't sufficiently far away from a body, depends on the quotient s/d where s is the width
of the region represented by the internal node, and d is the distance between the body and the node's center of
mass. The node is sufficiently far away when this ratio is smaller than a threshold value θ. The parameter θ
determines the accuracy of the simulation; larger values of θ increase the speed of the simulation but
decreases its accuracy. If θ = 0, no internal node is treated as a single body and the algorithm degenerates
to a direct-sum algorithm.


start at root node
if node is external, calculate force exerted by node on body
    if node is internal, calculate s/d
        if s/d < θ, treat internal node as single body, calculate force
        if s/d > θ, recursively apply Barnes-Hut to each child node
        add forces together
'''






# create some particles
particles = [
    Particle(
        pos=[random.uniform(-100, 100) for _ in range(3)],  # random position
        vel=[random.uniform(-1, 1) for _ in range(3)],  # random velocity
        mass=random.uniform(1, 1000),  # random mass
        radius=random.uniform(0.1, 1),  # random radius
        charge=random.uniform(-1, 1),  # random charge
        color=(random.random(), random.random(), random.random())  # random color
    ) 
    for _ in range(10000)
]

tree = BHTree(particles)
tree.root.treeview()












