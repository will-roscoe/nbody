from __future__ import annotations
import cProfile
import re

from numba import jit 

import gc
import math

import itertools


#import rich.traceback
#rich.traceback.install()


G = 6.67430e-11  # gravitational constant

import sys
from typing import Iterable, List, Union, Any
import numpy as np
import random

@jit(nopython=True,cache=True)
def dist(a: np.ndarray, b: np.ndarray) -> float:
    """calculate distance between two points"""
    return np.linalg.norm(a - b)


@jit(nopython=True)
def force_ab(pos1: np.ndarray, pos2: np.ndarray, mass1: np.float64, mass2: np.float64) -> np.ndarray:
    '''calculate force on target node due to trial node using Newton's law of gravitation'''
    f = np.clip(G*mass1*mass2/(dist(pos1,pos2)**3)*(pos1-pos2), -10**3, 10**3)
    return f

@jit(nopython=True)
def c_o_m(positions:np.ndarray, masses:np.ndarray) -> np.ndarray:
    total_mass = np.sum(masses, dtype=np.float64)
    total_pos = np.sum(positions*masses[:,None], axis=0, dtype=np.float64)
    return total_pos/total_mass

@jit(nopython=True)
def _contains(pos: np.ndarray, region: np.ndarray) -> bool:
    '''check if particle is within node's volume by comparing each dimension'''
    return np.all(np.logical_and(region[0] <= pos, pos <= region[1]))


##########################################################################################
#                                       Particle class                                   #
##########################################################################################


class Particle:
    def __init__(self,
                 pos: np.ndarray, vel: np.ndarray,
                 mass: float, radius: float,
                 charge: float, color: Union[str, np.ndarray] = 'white'):
        
        self.pos = np.array(pos, dtype=np.float64)
        self.vel = np.array(vel, dtype=np.float64)
        self.mass = mass
        self.radius = radius
        self.charge = charge
        self.color = color
        self.node = None
        self.p_id = 'None'
        
        self._pos = np.array(pos, dtype=float) # next position, where the mproc worker will put it.
        self._vel = np.array(vel, dtype=float)
        self.has_updates = False # flag for whether the particles _pos, _vel has been updated by the worker yet. assuming faster than checking if _pos == pos
        #debug print(f'Particle created at {np.round(self.pos,3)}')
        
    
    def update_changes(self,keeprefs=False) -> None:
        '''changes pending values to current values, and resets the flag''' 
        if self.has_updates:
            self.pos = self._pos
            self.vel = self._vel
            if not keeprefs:
                self.node.particles = []
                self.node = None
            self.has_updates = False
        else:
            raise ValueError("No updates to apply, but update_changes was called.")
    

    def __repr__(self) -> str:
        return f"Particle[{self.p_id}](position={np.round(self.pos,3)}, velocity={np.round(self.vel,3)},mass={np.round(self.mass)}), !UP:{int(self.has_updates)}, {f'dR:{int(np.all(self._pos == self.pos))}, dV:{int(np.all(self._vel == self.vel))}'}"
    

    def dist_p2p(self, other: Particle) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""
        return dist(self.pos, other.pos)
    
    def force_ab_p(self, trial: Union[Particle, Node]) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return force_ab(self.pos, trial.pos, self.mass, trial.mass)

    def collect_data(self) -> np.ndarray:
        return self.p_id,self.pos

##########################################################################################
#                                       Node class                                       #
##########################################################################################


class Node:
    iter_map = np.array(list(itertools.product([0, 1], repeat=3)))
    
    def __init__(self, region:np.ndarray[np.ndarray] , mass: float=0., depth=0, parent_tree: BHTree = None) -> None:
        self.particles = []  # particles directly in this node
        self.is_end = False # is this a leaf node
        self.depth = depth # depth of node in tree
        self.parent_tree = parent_tree
        self.children = []
        
        self.region = region # minimum and maximum coordinates of node, defining its volume.
        self.pos =  np.mean(self.region, axis=0, dtype=np.float64) # center of node
        self.lengths = self.region[1]-self.region[0] # each node is cubic so only need one size value
        self.size = np.linalg.norm(self.lengths) # maximum distance between two vertices of the node

        self.cmass = self.pos # center of mass of node (start at center of node for now)
        self.mass = mass # total mass of node, including children
        
        #debug print(f"{self.desc()}: created with depth {depth}:")
    
    def get_child_regions(self) -> Iterable[np.ndarray]:
        '''generate regions for all 8 children of this node'''
        childlens = (self.region[1]-self.region[0])/2 
        base_ =np.array((self.region[0],self.region[0]+childlens))
        for n in childlens*Node.iter_map:
            yield base_ + np.array([n,n])
    
    def insert(self, particle: Particle) -> None:
        '''insert particle into node's list of children and update mass of node'''
        self.particles.append(particle) 
        self.mass += particle.mass
        #debug print(f'{self.desc()}: particle {particle.p_id} inserted into node')

    def split(self) -> None:
        '''split node into 8 octants and redistribute particles to children nodes'''
        self.compute_cmass()
        if not self.is_end:
            # Create child nodes
            for region in self.get_child_regions():
                self.children.append(Node(region=region, depth=self.depth + 1, parent_tree=self.parent_tree))
            #debug print(f'{self.desc()}: all children created, beginning redistribution...')
            for particle in self.particles: # redistribute particles to children nodes
                #debug print(f'{self.desc()}: searching for node to insert particle {particle.p_id} at {np.round(particle.pos,4)}')
                for child in self.children:  
                    if child.contains(particle): # if particle is within child node, insert
                        child.insert(particle)
                        break        
                else:
                    raise ValueError(F'particle {particle.p_id} is not within any child node.')
                
            #debug print(f'{self.desc()}: particles redistributed, clearing lists')        
            self.children = [child for child in self.children if child.mass > 0] # get rid of empty nodes
            self.particles = []  # clear particles list after redistribution, since they are now in children
            
            for child in self.children:
                if len(child.particles) > 1: # if child has more than one particle, split again
                    child.split()
                else:
                    child.is_end = True # if child has only one particle, it is a leaf node
                    child.particles[0].node = child
                    child.cmass = child.particles[0].pos
                self.parent_tree.nodes.append(child)


    def contains(self, particle: Particle) -> bool:
        """check if particle is within node's volume"""
        return _contains(particle.pos, self.region)
    # check if particle is within node's volume by comparing each dimension 
    
    
    def compute_cmass(self) -> np.ndarray:
        """calculate center of mass of node using sum of mass*position / total mass. Does not need to be calculated after intitial particles have been inserted. use self.cmass instead. """
        positions, masses = np.array([p.pos for p in self.particles]), np.array([p.mass for p in self.particles])
        self.cmass = c_o_m(positions, masses)
    

    def dist_n2n(self, other: Node) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""   
        return dist(self.cmass, other.cmass)
    
    
    def force_ab_n(self, trial: Node) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return force_ab(self.cmass, trial.cmass, self.mass, trial.mass)

    
    # ----- Debugging methods -----
    def __repr__(self) -> str: # a debug representation of the node containing its children's and particles' info recursively
        return f'{self.desc()}: lengths={np.round(self.lengths,3)}, mass={np.round(self.mass)}, region={np.round(self.region,3).tolist()},'
    def __str__(self) -> str: # human readable representation of the node
        return (f"""{self.desc()}: region={np.round(self.region, 3).tolist()} lengths={np.round(self.lengths,3)} mass={round(self.mass,4)}, Leaf={int(self.is_end)}""")
    def desc(self):
        return f'd{self.depth}NODE[{len(self.particles)}p,{len(self.children)}ch]'
    def treeview(self) -> str: # recursive representation of the tree
        """print tree structure of node and children recursively"""
        bars = "   ⎸"*self.depth
        info = f"⟶  d{self.depth} Node[{len(self.children)}ch,{len(self.particles)}p], Region:{np.round(self.region,3)}, Mass:{round(self.mass,4)}"
        if not self.is_end:
            info = "\u001b[37m\u001b[1m" + info + "\u001b[37m\u001b[0m"
        print(bars + info)
        for c in self.children:
            c.treeview()
    


##########################################################################################
#                                        BHTree class                                    #
##########################################################################################

class BHTree:
    _off = np.array([0.001,]*3) # adding a small offset to the region to prevent particles from being on the edge of the region, which causes errors in the tree. completely arbitrary value.
    def __init__(self, particles=None) -> None:
        self.nodes = self.update(particles)
        self.root = self.nodes[0]


    def _cull(self) -> None:
        '''clears all nodes references'''
        del self.nodes, self.root
        
        
    def get_new_root(self, particles: List[Particle]) -> Node:
        positions = np.array([p.pos for p in particles])
        new_region = np.array((np.min(positions, axis=0)-BHTree._off,np.max(positions, axis=0)+BHTree._off))
        new_root = Node(region=new_region, parent_tree=self) 
        return new_root
    
    def insert_particles(self, particles: List[Particle], root_node: Node) -> Node:
        for particle in particles:
            if not root_node.contains(particle):
               raise ValueError("particle is not within root node region before tree initialization.")
            root_node.insert(particle) # insert particles into root node 
        return root_node
    
    
    def update(self, particles: List[Particle]) -> List[Node]:
        self.root = self.insert_particles(particles, self.get_new_root(particles))
        self.nodes = [self.root] # nodes list contains list of nodes
        self.root.split() # begin recursive splitting of nodes
        return self.nodes
    
##########################################################################################
#                                      ParticleArtist class                              #
########################################################################################## 

class DataCollector:
    def __init__(self, engine: Engine) -> None:
        self.engine = engine
        self.particles = engine.particles
        self.data = {p.p_id:np.array(p.pos, ndmin=2) for p in self.particles}
    
    def collect_data(self) -> None:
        next_data = [p.collect_data() for p in self.particles]
        for i,pos in next_data:
            self.data[i] = np.vstack((self.data[i],pos))



  
##########################################################################################
#                                      Engine class                                      #
########################################################################################## 

class Engine:
    def __init__(self, particles: List[Particle], dt: float = 1e-3, theta: float=0.5, recorder=DataCollector) -> None:
        #debug print('ENG: init start')
        for i, particle in enumerate(particles):
            particle.p_id = i
        self.particles = particles
        self.dt = dt
        self.tree = BHTree(particles)
        self.theta = theta
        self._data = recorder(self)
        
    def update(self) -> None:
        self._data.collect_data()
        for particle in self.particles:
            particle._pos += particle.vel * self.dt
            particle._vel +=  self.force_on(particle) * self.dt / particle.mass
            particle.has_updates = True

        self.tree._cull()
        for particle in self.particles:
            particle.update_changes()
        self.tree.update(self.particles)

    
    #?####################################################################################
    #?                               Barnes-Hut algorithm part                          ##
    #?####################################################################################
    #@jit(forceobj=True)
    def force_on(self, target_particle:Particle) -> np.ndarray:
        '''get the force on a target particle from the tree using the Barnes-Hut algorithm'''
        target = target_particle.node
        trials = [self.tree.root] 
        candidates = []
        while trials:
            trial = trials.pop(0)
            if target != trial: # if target is trial, then throw away
                if trial.is_end: # if it is leaf node, calculate force
                    candidates.append(trial)
                else:
                    if  trial.size/trial.dist_n2n(target) < self.theta: # if condition is met for branch node, add to list
                        candidates.append(trial)
                    else: # otherwise add children to trials
                        trials.extend(trial.children)
            else: # throw away if none of the above conditions are met
                pass       
        return np.sum([target.force_ab_n(trial) for trial in candidates], axis=0)
    



































##########################################################################################
#                                           Main                                         #
########################################################################################## 

def main():
    # create some particles

    print('initializing particles')
    random.seed(10000)
    particles_ = [
        Particle(
            pos=[random.uniform(-100, 100) for _ in range(3)],  # random position
            vel=[random.uniform(-100, 100) for _ in range(3)],  # random velocity
            mass=random.uniform(30, 1000),  # random mass
            radius=random.uniform(0.1, 1),  # random radius
            charge=random.uniform(-1, 1),  # random charge
            color=(random.random(), random.random(), random.random())  # random color
        ) 
        for _ in range(10)
    ]
    
    if __name__ == '__main__':
        particles = [Particle(pos=[1,0,0], vel=[1,0,0], mass=100, radius=1, charge=0, color='red'),
                     Particle(pos=[-1,0,0], vel=[-1,0,0], mass=100, radius=1, charge=0, color='blue')]
        eng = Engine(particles, dt=0.1)
        for _ in range(1000):
            eng.update()
            print(_)
        print(eng._data.data)
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        plt.style.use('dark_background')
        for d in eng._data.data.values():
            ax.plot(*d.T)
        plt.show()
main()




