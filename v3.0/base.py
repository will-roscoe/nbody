
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
##########################################################################################
#                                       Particle class                                   #
##########################################################################################


class Particle:
    def __init__(self,
                 pos: np.ndarray, vel: np.ndarray,
                 mass: float, radius: float,
                 charge: float, color: Union[str, np.ndarray] = 'white'):
        
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

##########################################################################################
#                                       Node class                                       #
##########################################################################################


class Node:
    def __init__(self, pos: np.ndarray, size: float, mass: float=0., depth=0) -> None:
        self.pos = pos # center of node
        self.size = size # each node is cubic so only need one size value
        self.mass = mass # total mass of node, including children
        self.children = []  # child nodes
        self.particles = []  # particles directly in this node
        self.is_end = False # is this a leaf node
        self.depth = depth # depth of node in tree
        self.cmass = pos # center of mass of node (start at center of node for now)
        #debug print(f"d{depth} Node: pos={self.pos}, size={self.size}, mass={self.mass}")
    
    
    def insert(self, particle: Particle) -> None:
        '''insert particle into node and update mass of node'''
        self.particles.append(particle) 
        self.mass += particle.mass
        #debug print(f"inserted particle with mass {particle.mass} into node with mass {self.mass}")

    
    def split(self) -> None:
        '''split node into 8 octants and redistribute particles to children nodes'''
        if not self.is_end: # if it is already a leaf node, don't split
            #debug print(f"splitting a d{self.depth} Node: pos={self.pos}, size={self.size}, mass={self.mass}")
            #split into 8 octants
            for i in [1,-1]:  # bisecting in each dimension
                for j in [1,-1]:
                    for k in [1,-1]:
                        pos = self.pos + np.array([i, j, k]) * self.size / 2 # calculate child node position
                        self.children.append(Node(pos, self.size/2, depth=self.depth+1))  # add as child node
            
            for particle in self.particles: # redistribute particles to children nodes
                for child in self.children:  
                    if child.contains(particle): # if particle is within child node, insert
                        child.insert(particle)
                        break
                else:
                    raise ValueError(f"particle not in any child node: \n{particle}")
                    
            self.children = [child for child in self.children if child.mass > 0] # get rid of empty nodes
            self.particles = []  # clear particles list after redistribution, since they are now in children
            for child in self.children:
                if len(child.particles) > 1: # if child has more than one particle, split again
                    child.split()
                    self.fast_cmass() # update center of mass of parent node
                else:
                    child.is_end = True # if child has only one particle, it is a leaf node
    
    
    def contains(self, particle: Particle) -> bool:
        """check if particle is within node's volume"""
        return np.all(np.abs(particle.pos - self.pos) <= self.size)  # check if particle is within node's volume by comparing each dimension 
    
    
    def fast_cmass(self) -> np.ndarray:
        """calculate center of mass of node using fast method"""
        self.cmass = np.sum([p.mass*p.pos for p in [*self.particles,*self.children]], axis=0) / self.mass # sum of mass*position / total mass. faster because we remove most children before using it.
        return self.cmass
    
    
    def __repr__(self) -> str: # a debug representation of the node containing its children's and particles' info recursively
        return (f"""
d{self.depth}Node: pos={self.pos} size={self.size} mass={self.mass} 
children({len(self.children)}): {[child for child in self.children]}
particles({len(self.particles)}): {[particle for particle in self.particles]}""")


    

    def treeview(self) -> str: # recursive representation of the tree
        """print tree structure of node and children recursively"""
        bars = "   ⎸"*self.depth
        info = f"⟶  d{self.depth} Node[{len(self.children)}ch,{len(self.particles)}p], cmass {self.cmass}, Mass:{self.mass}"
        if not self.is_end:
            info = "\u001b[37m\u001b[1m" + info + "\u001b[37m\u001b[0m"
        print(bars + info)
        for c in self.children:
            c.treeview()

    
    def dist(self, other) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""
        return np.linalg.norm(self.cmass - other.cmass)

##########################################################################################
#                                        BHTree class                                    #
##########################################################################################

class BHTree:
    def __init__(self, particles, theta:float=0.5, origin:np.ndarray=np.zeros(3), size:float=100.) -> None:
        self.theta = theta
        #debug [np.mean(comp) for comp in zip(*[p.pos for p in particles])]
        #debug initpos = np.array([np.mean(comp) for comp in zip(*[p.pos for p in particles])])
        self.root = Node(
            pos= origin, size = size    #commented out bits are for a case where we might not know the size of the system initially
            #debug pos = initpos,
            #debug size = np.max([np.max([np.abs(particle.pos - initpos)]) for particle in particles])
            )
        for particle in particles:
            self.root.insert(particle) # insert particles into root node 
        self.root.split() # begin recursive splitting of nodes


        
        

    #?####################################################################################    

##########################################################################################
#                                      Engine class                                      #
########################################################################################## 




##########################################################################################
#                                           Main                                         #
########################################################################################## 


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












