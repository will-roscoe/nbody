from __future__ import annotations
import cProfile
from calendar import c
import math

import itertools


# Same as add_hook(always=True)
import colored_traceback.always





G = 6.67430e-11  # gravitational constant

import sys
from typing import Iterable, List, Union, Any
import numpy as np
import random

from vispy import scene
from vispy.scene import visuals
from vispy.color import Color


import multiprocessing as mproc

from queue import Queue

from pympler.asizeof import asizeof

from time import perf_counter as perftime
from time import time
import gc
import concurrent.futures 
update_time = perftime()
def ptime(start=update_time): return round(perftime()-start, 5)

DEBUG_PRINTING = True
LOG_MSG = True
DISPLAY_TIME = True
logging_cache = []
log_dir = 'logs'
log_file = f'{log_dir}/log{int(time())}.txt'
open(log_file, 'w').close()

def msg(message: 'str', obj: Any = None) -> None: 
    time = (f'{round(ptime(),4)}s:'if DISPLAY_TIME else '')
    fmtmsg = f'{str(time) : >8}{" "+str(id(obj)): <20}{str(message): <100}'
    if DEBUG_PRINTING: print(fmtmsg)
    if LOG_MSG:
        with open(log_file, 'a') as f:
            f.write(fmtmsg+'\n')

def err(errtype=ValueError, message: str = 'An Error Has occurred', reporter: Any = None, target: Any = None) -> None:
    if reporter is None:
        reporter = 'Not Specified'
    if target is None:
        target = 'Not Specified'
    msg(f'{str(errtype)}: {message}, Rep:{repr(reporter)}, Tar:{repr(target)}', reporter)
    raise errtype(f'{message}\nRaised by: {repr(reporter)}\nTargeting: {repr(target)}')
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
        self.node = None
        self.p_id = 'None'
        
        self._pos = np.array(pos, dtype=float) # next position, where the mproc worker will put it.
        self._vel = np.array(vel, dtype=float)
        self.has_updates = False # flag for whether the particles _pos, _vel has been updated by the worker yet. assuming faster than checking if _pos == pos
    def check_pos_vel_changes(self) -> None:
        '''checks if there are changes to the position and velocity of the particle'''
        return f'P=_P:{np.all(self._pos == self.pos)}, V=_V:{np.all(self._vel == self.vel)}'
        
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
            err(ValueError, f"No updates to apply, but update_changes was called.", self, self)
    
    def __repr__(self) -> str:
        return f"Particle[{self.p_id}](position={self.pos}, velocity={self.vel},mass={self.mass}), has_updates={self.has_updates}, {self.check_pos_vel_changes()}"
    
    def dist_p2p(self, other: Particle) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""
        return np.linalg.norm(self.pos - other.pos)
    
    def force_ab(self, trial: Union[Particle, Node]) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return G*self.mass*trial.mass/(self.dist(trial)**3)*(self.pos-trial.pos)




##########################################################################################
#                                       Node class                                       #
##########################################################################################


class Node:
    def __init__(self, region:np.ndarray[np.ndarray] , mass: float=0., depth=0, parent_tree: BHTree = None) -> None:
        self.particles = []  # particles directly in this node
        self.is_end = False # is this a leaf node
        self.depth = depth # depth of node in tree
        self.parent_tree = parent_tree
        self.children = []
        
        self.region = region # minimum and maximum coordinates of node, defining its volume.
        self.pos =  np.mean(self.region, axis=0) # center of node
        self.lengths = self.region[1]-self.region[0] # each node is cubic so only need one size value
        self.cmass = self.pos # center of mass of node (start at center of node for now)
        self.mass = mass # total mass of node, including children
        msg(f"{self.desc()}: created with depth {depth}: {str(self)}")
    
    
    def insert(self, particle: Particle) -> None:
        '''insert particle into node and update mass of node'''
        
        msg(f'{self.desc()}: inserting particle {particle.p_id} into this node',self)
        self.particles.append(particle) 
        self.mass += particle.mass
        msg(f'{self.desc()}: particle {particle.p_id} inserted into node',self)

    
    def split(self) -> None:
        '''split node into 8 octants and redistribute particles to children nodes'''
        msg(f'{self.desc()}: node split start',self)
        self.fast_cmass()
        
        if not self.is_end:
            child_lengths = self.lengths/2 # if it is already a leaf node, don't split
            msg(f'{self.desc()}: splitting node...',self)
            #debug print(f"splitting a d{self.depth} Node: pos={self.pos}, size={self.size}, mass={self.mass}")
            #split into 8 octants
            child_lengths = self.lengths / 2
            half_child_lengths = child_lengths / 2
            child_positions = self.pos + half_child_lengths * np.array(list(itertools.product([-1, 1], repeat=3)))
            children_regions = np.array([child_positions + offset * child_lengths for offset in itertools.product([0, 1], repeat=3)])
            # Create child nodes
            for region in children_regions:
                self.children.append(Node(region=region, depth=self.depth + 1, parent_tree=self.parent_tree))
           
            msg(f'{self.desc()}: all children created, beginning redistribution...',self)
            for particle in self.particles: # redistribute particles to children nodes
                
                msg(f'{self.desc()}: searching for node to insert particle {particle.p_id} at {np.round(particle.pos,4)}',self)
                for child in self.children:  
                    msg(f'{self.desc()}: testing child node {str(child)}',self)
                    
                    if child.contains(particle): # if particle is within child node, insert
                        msg(f'{self.desc()}: child node found. inserting particle {particle.p_id} into child node {str(child)}',self)
                        child.insert(particle)
                        msg(f'{self.desc()}: inserted particle {particle.p_id} into node',self)
                        break        
                else:
                    err(ValueError, f'{particle.p_id} is not within any child node.', self, particle)
                
            msg(f'{self.desc()}: particles redistributed, clearing lists',self)        
            self.children = [child for child in self.children if child.mass > 0] # get rid of empty nodes
            self.particles = []  # clear particles list after redistribution, since they are now in children
            msg(f'{self.desc()}: lists cleared, starting child splits',self)
            
            for child in self.children:
                msg(f'{self.desc()}: querying child node {child.desc()} or further splitting: ',self)
                 # update center of mass of this node
                if len(child.particles) > 1: # if child has more than one particle, split again
                    msg(f'{self.desc()}: child node {child.desc()} has more than one particle, asking it to split',self)    
                    child.split()
                    msg(f'{self.desc()}: child node {str(child)} finished splitting',self)
                else:
                    msg(f'{self.desc()}: node {str(child)} is a leaf, updating its attributes',self)
                    child.is_end = True # if child has only one particle, it is a leaf node
                    child.particles[0].node = child
                    child.cmass = child.particles[0].pos
                    msg(f'{self.desc()}: node {str(child)} is updated',self)
                msg(f'{self.desc()}: adding child node {str(child)} to tree\'s list',self)
                self.parent_tree.nodes.append(child)

            
            msg(f'{self.desc()}: node split end',self)
        

    def contains(self, particle: Particle) -> bool:
        """check if particle is within node's volume"""
        msg(f'{self.desc()}: checking if particle {particle.p_id} is within node',self)
        return np.all(np.logical_and(self.region[0] <= particle.pos, particle.pos <= self.region[1]))
    # check if particle is within node's volume by comparing each dimension 
    
    
    def fast_cmass(self) -> np.ndarray:
        """calculate center of mass of node using fast method"""
        #debug msg(f'{self.desc()}: calculating center of mass of node',self)
        self.cmass = np.sum([*[p.mass*p.pos for p in self.particles],*[p.cmass*p.mass for p in self.children]], axis=0) / self.mass # sum of mass*position / total mass. faster because we remove most children before using it.
        return self.cmass
    
    
    def __repr__(self) -> str: # a debug representation of the node containing its children's and particles' info recursively
        return (f"""d{self.depth}Node: region={self.region}, lengths={self.lengths}, mass={self.mass}, 
children({len(self.children)}): \n {[child for child in self.children]}
particles({len(self.particles)}): \n {[particle for particle in self.particles]}""")


    def __str__(self) -> str: # human readable representation of the node
        return (f"""d{self.depth}Node: region={np.round(self.region, 3)} lengths={np.round(self.lengths,3)} mass={round(self.mass,4)} [{len(self.children)}ch,{len(self.particles)}p], Leaf={self.is_end}""")
    
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
    
    
    def dist(self, other: Node) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""
        #debug msg(f'{self.desc()}: calculating distance between this node and {str(other)}',self)   
        return np.linalg.vector_norm(self.cmass - other.cmass)
    
    def force_ab(self, trial: Node) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return G*self.mass*trial.mass/(self.dist(trial)**3)*(self.cmass-trial.cmass)


##########################################################################################
#                                        BHTree class                                    #
##########################################################################################

class BHTree:
    def __init__(self, particles=None) -> None:
        msg('BHT: init start',self)
        self.nodes = self.update(particles)
        msg('BHT: init',self)
        self.root = self.nodes[0]


    def begin_cull(self) -> None:
        '''clears all particles from nodes and sets mass to 0'''
        msg('BHT: culling nodes start',self)
        msg(f'BHT: currently referenced by:{[id(n) for n in gc.get_referrers(*self.nodes)]}')
        while self.nodes:
            node = self.nodes.pop(0)
            node.mass = 0
            node.particles = []
            node.children = []
            del node
        msg(f'BHT: currently referenced by:{[id(n) for n in gc.get_referrers(*self.nodes)]}')
        msg('BHT: culled',self)
        
        
    def get_new_root(self, particles) -> Node:
        #pos_ = np.array([np.mean(comp) for comp in zip(*[p.pos for p in particles])])
        positions = np.array([p.pos for p in particles])
        new_region = np.array((np.min(positions, axis=0),np.max(positions, axis=0)))
        new_root = Node(region=new_region, parent_tree=self) 
        msg(f'BHT: root node created: {str(new_root)}, sending to update',self)
        return new_root, new_region
    
    def update(self, particles) -> None:
        msg('BHT: tree update start',self)
        self.root, root_region = self.get_new_root(particles)
        for particle in particles:
            msg(f'BHT: inserting particle [{particle.p_id}] into root node',self)
            if not self.root.contains(particle):
                err(ValueError, f"particle is not within root node region before tree initialization.", self, particle)
            self.root.insert(particle) # insert particles into root node 
        self.nodes = [self.root] # nodes list contains list of nodes
        self.root.split() # begin recursive splitting of nodes
        msg(f"BHT: tree updated with {len(self.nodes)} total nodes",self)
        msg(f"{self.root.treeview()}...",self)
        return self.nodes
    def vispy_draw(self) -> None:
        canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)

        # Set up a viewbox to display the cube with interactive arcball
        view = canvas.central_widget.add_view()
        view.bgcolor = '#efefef'
        view.camera = 'turntable'
        view.padding = 100

        color = Color("#3f51b5")
        queued_nodes = [self.root]
        while queued_nodes:
            node = queued_nodes.pop(0)
            cube = scene.visuals.Box(planes = (0,1,0,1,0,1), color=color, edge_color="black",
                                parent=view.scene)
            if not node.is_end:
                queued_nodes.extend(node.children) 
        
        if __name__ == '__main__' and sys.flags.interactive == 0:
            canvas.app.run()

    
##########################################################################################
#                                      Engine class                                      #
########################################################################################## 

class Engine:
    def __init__(self, particles: List[Particle], dt: float = 1e-3, theta: float=0.5) -> None:
        msg('ENG: init start',self)
        for i, particle in enumerate(particles):
            particle.p_id = i
        self.particles = particles
        self.dt = dt
        self.tree = BHTree(particles)
        self.theta = theta
        msg('ENG: init',self)
        
    def update(self) -> None:
        msg('ENG: update start',self)
        
        for particle in self.particles:
            msg(f'ENG: updating particle {particle.p_id}',self)
            particle._pos += particle.vel * self.dt
            particle._vel +=  self.force_on(particle) * self.dt / particle.mass
            particle.has_updates = True
        self.tree.begin_cull()
        for particle in self.particles:
            msg(f'ENG: updating live values of {particle.p_id}',self)
            particle.update_changes()
        self.tree.update(self.particles)
        msg('ENG: update success',self)
    
    #?####################################################################################
    #?                               Barnes-Hut algorithm part                          ##
    #?####################################################################################
    
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
                    if  trial.size/trial.dist(target) < self.theta: # if condition is met for branch node, add to list
                        candidates.append(trial)
                    else: # otherwise add children to trials
                        trials.extend(trial.children)
            else: # throw away if none of the above conditions are met
                pass       
        return np.sum([target.force_ab(trial) for trial in candidates], axis=0)


##########################################################################################
#                                           Main                                         #
########################################################################################## 


# create some particles
msg('initializing particles')
random.seed(10000)
particles_ = [
    Particle(
        pos=[random.uniform(-100, 100) for _ in range(3)],  # random position
        vel=[random.uniform(-1, 1) for _ in range(3)],  # random velocity
        mass=random.uniform(1, 1000),  # random mass
        radius=random.uniform(0.1, 1),  # random radius
        charge=random.uniform(-1, 1),  # random charge
        color=(random.random(), random.random(), random.random())  # random color
    ) 
    for _ in range(10)
]
msg('particles initialized')





if __name__ == '__main__':
    particles = particles_[:]
    msg('initializing engine')
    eng = Engine(particles)
    msg('starting update loop')
    p1 = ptime()
    def pitime(start=p1): print(f'T={round(ptime()-start, 5)}s')
    for _ in range(100):
        eng.update()
        pitime()



    
    
    
def test():   
    runs = []
    for y in range(1,11):
        print(f"running test with {y} particles")
        for x in range(3):
            particles = particles_[:5000]
            #print(asizeof(particles))
            p1c = ptime()
            engine = Engine(particles)
            p1a = ptime()
            [engine.tree.force_on(particles[m]) for m in range(y)]
            p1b = ptime()
            #print(f'N#{x} ENG~{asizeof(engine)}bytes')
            
            particles = particles_[5000:]
            #print(asizeof(particles))
            p2c = ptime()
            engine = Engine(particles, tree=False)
            p2a = ptime()
            [engine.force_on_old(particles[n]) for n in range(y)]
            p2b = ptime()
            #print(f'O#{x} ENG~{asizeof(engine)}bytes')
            runs.append([p1c, p1a, p1b, p2c, p2a, p2b, x,y])
        #for run in runs: 
            #print(run)
        r1, r2, r3, r4 = [sum([run[i] for run in runs])/len(runs) for i in range(4)]
        print(f"avgs: (Engine&/Tree Creation + Calc) \nN:{round(r3, 5)+round(r1, 5)}s \nO:{round(r4, 5)+round(r2, 5)}s")
        print(f"speedup on calc: {r1/r2} (N/O) \nspeedup on init: {r3/r4} (N/O) \nspeedup overall: {(r1+r3)/(r2+r4)}\n(bigger is worse)")
        
        import csv
        with open('output.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(runs)





