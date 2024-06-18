from __future__ import annotations
# numpy based vector objects, fixed dtype float64

# record vectors need to have variable length, but numpy won't allow this with np.array .
# we could use a list of numpy arrays but it may be slow to access elements.
# we need to be able to acess elements quickly within a numpy system with variable length, so 
G = 6.67430e-11  # gravitational constant
k = 8.9875517923e9  # Coulomb's constant  # permeability of free space


import sys
from typing import List, Union
import numpy as np
import random

from vispy import scene
from vispy.scene import visuals
from vispy.color import Color

import tqdm as tq

import multiprocessing as mproc
mproc.set_start_method('spawn', force=True)

from pympler.asizeof import asizeof

from time import perf_counter as perftime

update_time = perftime()
def ptime(start=update_time): return round(perftime()-start, 5)




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
    
    def update_changes(self):
        if self.has_updates:
            self.pos = self._pos
            self.vel = self._vel
            self.has_updates = False
        else:
            raise ValueError(f"particle {self.p_id} has no updates to apply, but update_changes was called.")
    def __repr__(self) -> str:
        return f"Particle(p={self.pos},v={self.vel},m={self.mass},r={self.radius},c={self.charge})\n"
    def dist(self, other) -> np.ndarray: 
        """calculate distance between center of masses of two nodes"""
        return np.linalg.norm(self.pos - other.pos)
    def force_ab(self, trial) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return G*self.mass*trial.mass/(self.dist(trial)**3)*(self.pos-trial.pos)





##########################################################################################
#                                       Node class                                       #
##########################################################################################


class Node:
    def __init__(self, pos: np.ndarray, size: float, mass: float=0., depth=0, parent_tree: BHTree = None) -> None:
        self.pos = pos # center of node
        self.size = size # each node is cubic so only need one size value
        self.mass = mass # total mass of node, including children
        self.children = []  # child nodes
        self.particles = []  # particles directly in this node
        self.is_end = False # is this a leaf node
        self.depth = depth # depth of node in tree
        self.cmass = pos # center of mass of node (start at center of node for now)
        self.parent_tree = parent_tree
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
                        self.children.append(Node(pos, self.size/2, depth=self.depth+1, parent_tree=self.parent_tree))  # add as child node
            
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
                else:
                    child.is_end = True # if child has only one particle, it is a leaf node
                    child.particles[0].node = child
                self.parent_tree.nodes.append(child)
            self.fast_cmass() # update center of mass of this node
        

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


    def __str__(self) -> str: # human readable representation of the node
        return (f"""d{self.depth}Node: pos={self.pos} size={self.size} mass={self.mass} [{len(self.children)}ch,{len(self.particles)}p]""")
    

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
    
    def force_ab(self, trial) -> np.ndarray:
        '''calculate force on target node due to trial node using Newton's law of gravitation'''
        return G*self.mass*trial.mass/(self.dist(trial)**3)*(self.cmass-trial.cmass)





##########################################################################################
#                                        BHTree class                                    #
##########################################################################################

class BHTree:
    def __init__(self, particles, theta:float=0.5, origin:np.ndarray=np.zeros(3), size:float=100.) -> None:
        self.theta = theta
        #debug [np.mean(comp) for comp in zip(*[p.pos for p in particles])]
        #debug initpos = np.array([np.mean(comp) for comp in zip(*[p.pos for p in particles])])
        self.root = Node(
            pos= origin, size = size,    #commented out bits are for a case where we might not know the size of the system initially
            #debug pos = initpos,
            #debug size = np.max([np.max([np.abs(particle.pos - initpos)]) for particle in particles])
            parent_tree=self)
        for particle in particles:
            self.root.insert(particle) # insert particles into root node 
        self.nodes = [self.root] # nodes list contains list of nodes at each depth level 
        self.root.split() # begin recursive splitting of nodes
        #print(f"tree initialized with {len(self.nodes)} total nodes")
    
    
    def vispy_draw(self):
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


    
    
    #?####################################################################################
    #?                               Barnes-Hut algorithm part                          ##
    #?####################################################################################
    def choose_nodes(self, target: Node, trial: Node) -> int:
        """chooses whether to:
        [0]:add children to trials or
        [1]:calculate force AB or
        [2]:do nothing"""
        if target == trial: # if target is trial, then throw away
            return 2
        if trial.is_end: # if it is leaf node, calculate force
            return 1
        elif not trial.is_end:
            if  trial.size/trial.dist(target) < self.theta: # if condition is met for branch node, calculate force
                return 1 
            else: # otherwise add children to trials
                return 0
        else: # throw away if none of the above conditions are met
            return 2
    

    def force_on(self, target_particle:Particle) ->np.ndarray:
        target = target_particle.node
            # picking relevant nodes by recursively checking conditions of self.choose_nodes on each node
        trials = [self.root]
        candidates = []
        while trials:
            trial = trials.pop(0)
            action = self.choose_nodes(target, trial)
            if action == 0:
                trials.extend(trial.children)
            elif action == 1:
                candidates.append(trial)
            elif action != 2:
                raise ValueError(f"invalid action: {action}")
        
        
        ''' may not be worth it to parallelize this part, had problems with is returning as nonetype parralelized again.
            might work with an approach using a pipe between processes but might not be worth it here
            https://shorturl.at/9bGZt Pipes Docs

        if __name__ == '__main__':
            pool = mproc.Pool(4)
            results = pool.starmap(self.force_ab, [(target, trial) for trial in candidates])
            pool.close()
            pool.join()
            return np.sum(results, axis=0)
        '''
        return np.sum([target.force_ab(trial) for trial in candidates], axis=0) #quick and dirty non-parallelized version
    

    #?####################################################################################    

##########################################################################################
#                                      Engine class                                      #
########################################################################################## 

class Engine:
    def __init__(self, particles: List[Particle], dt: float = 1e-3, tree=True):
        for i, particle in enumerate(particles):
            particle.p_id = i
        self.particles = particles
        self.dt = dt
        self.tree = [BHTree(particles) if tree else None][0]

    def force_on_old(self, target_particle:Particle) ->np.ndarray:
        candidates = [p for p in self.particles if p != target_particle]     
        return np.sum([target_particle.force_ab(trial) for trial in candidates], axis=0)

    def update_old(self) -> None:
        for particle in self.particles:
            particle._pos = particle.pos + particle.vel * self.dt
            particle._vel = particle.vel + self.force_on_old(particle) * self.dt / particle.mass
            particle.has_updates = True
            #print(f'particle {particle.p_id} updated')
        for particle in self.particles:
            particle.update_changes()


    def update(self) -> None:
        for particle in self.particles:
            particle._pos = particle.pos + particle.vel * self.dt
            particle._vel = particle.vel + self.tree.force_on(particle) * self.dt / particle.mass
            particle.has_updates = True
            #print(f'particle {particle.p_id} updated')
        for particle in self.particles:
            particle.update_changes()

    def update_particle(self, particle: Particle, tree: BHTree) -> None:
        particle._pos = particle.pos + particle.vel * self.dt
        particle._vel = particle.vel + tree.force_on(particle) * self.dt / particle.mass
        particle.has_updates = True
        #print(f'particle {particle.p_id} ready to be updated')
    
    
    def pool_update(self) -> None:
        pool = mproc.Pool(4)
        results = pool.starmap(self.update_particle, [(particle, self.tree) for particle in self.particles])
        pool.close()
        pool.join()
        for particle in self.particles:
            particle.update_changes()
    '''
    def updater_worker(self, old_particles:mproc.Queue, new_particles:mproc.Queue) -> None:
        pid = mproc.current_process().pid
        new_values = []
        print(f'worker {pid} started, waiting for new particles...')
        while not old_particles.empty(): 
            print(f'{pid}: awaiting particle')
            particle = old_particles.get(timeout=1)
            p_id = particle.p_id
            print(f'{pid}: got particle {p_id}')
            _pos = particle.pos + particle.vel * self.dt
            _vel = particle.vel + self.tree.force_on(particle) * self.dt / particle.mass
            has_updates = True
            new_values.append((particle,p_id, _pos, _vel, has_updates))
            print(f'{pid}: finished with particle {p_id}')
            #new_particles.put(particle, timeout=1)
            #print(f'{pid}: returned particle {p_id}')
            #job_count += 1
        old_particles.close()
        #new_particles.close()
        job_count = len(new_values)
        print(f'{pid}: finished updating {job_count} particles, closed, exiting...')

    def update(self) -> None:
        awaiting_update = mproc.Queue()
        finished_update = mproc.Queue()
        workers = []
        for particle in self.particles:
            awaiting_update.put(particle)
        
        for _ in range(1):#mproc.cpu_count()//2):
            proc = mproc.Process(target=self.updater_worker, args=(awaiting_update, finished_update))
            workers.append(proc)
        print(f'MAIN: created {len(workers)} worker instances')
        for proc in workers:
            proc.start()
            print(f'MAIN: started worker {proc.pid}')
        print(f'MAIN: {len(mproc.active_children())} worker instances are active.')
        
        for proc in mproc.active_children():
            print(f'MAIN: waiting for worker {proc.pid}')
            proc.join()
            print(f'MAIN: worker {proc.pid} joined')
        print(f'MAIN: all workers joined, replacing self.particles')
        for _ in range(len(self.particles)):
            p = finished_update.get()
            p.update_changes()
        self.tree = BHTree(self.particles)
        print(f'MAIN: updated particles and tree')
    '''

##########################################################################################
#                                           Main                                         #
########################################################################################## 


# create some particles
random.seed(0)
particles_ = [
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






if __name__ == '__main__':
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







