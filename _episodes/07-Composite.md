---
title: 'Composite Design Pattern'
teaching: 30
exercises: 0
questions:
- 'How can an object notify its own state to an open-ended number of objects?'
keypoints:
- 'The composite design pattern provides a way for manipulating a tree data structure'
objectives:
- 'Show an example of the composite design pattern relavant to the Computational Molecular Sciences domain.' 
- 'Understand the composite design pattern.'
questions:
- 'How to compose objects into tree structures?'
- 'How to define a common interface for branches and leaves within a tree?'
---


## Problem

Consider a Monte Carlo (MC) simulation of a strongly
associating system. Conducting a 'naive' 
simulation on such systems using conventional
MC algorithms would require prohibitively amounts of time, as
the aggregate structures do not change sufficiently often
during the course of a simulation (i.e. clusters are not created or destroyed
frequently). To reproduce the correct statistics, it is necessary to 
design and use algorithms that more efficiently sample
configurations of the ensemble of interest. Cluster moves are one tool
to give more rapid sampling of aggregate configurations. 

In a typical cluster move, you might find two kinds of objects: particles and 
clusters. Particles are elements that tend to aggregate
into clusters. For instance, a surfactant molecule that tend to aggregate into larger
micellar structures. A cluster object usually contains one or more particles 
and thus is a composite
object, while a particle is a simple object. Cluster and particles may have similar
operations and attributes, such as center-of-mass translation or solid body rotation. It
would be convenient treat both cluster and
particle objects uniformly by defining a common interface.

A basic
cluster algorithm would involve the following steps

1. Check for the formation of clusters
2. Select a cluster randomly
3. Make a random perturbation using the selected cluster, such as
translation, rotation, removal or
addition of particles to the selected cluster.
4. Accept or reject the new configuration.

In this example, we will imagine we are implementing functionality to address point 3. 
We would like
to create objects for particles and clusters and provide them with methods for 
center-of-mass translation. 

## Solution: the composite design pattern

The composite design pattern allows you to create tree structures
and handle their components (branches and leaves) uniformly.
The composite pattern uses the following elements

1. Component. This is an abstract class that will serve
as parent for the branch and leaf elements. In our example, a
Component class would be an abstraction for patchy particles (leaves)
and clusters (branches).

~~~
from abc import ABC, abstractmethod

class Component(ABC):

    '''Base class for a Component. This is the 'component'
    base class in the pattern.'''

    @abstractmethod
    def translate(self):
        pass

~~~
{: .language-python}

2. Leaf. These are objects that have no children. They implement methods 
described by the Component parent class. Patchy particles, LJ particles
or Gay Berne particles might serve as leaf objects.

~~~
class Patchy(Component):

    '''A patchy particle is a type of particle.'''

    def translate(self):
        print('translating patchy particle')
~~~
{: .language-python}

3. Branch. This element stores child components and implement methods
defined by the component interface.

~~~
class Cluster(Component):

    '''A cluster is a collection of particles. It can
    be composed of spherical, anisotropic, patchy or 
    any other particle. This is a 'composite' in the
    Composite patern'''

    def __init__(self):
        self._particles = []

    def translate(self):
        for children in self._particles:
            children.translate()

    def add_particle(self, particle):
        self._particles.append(particle)
~~~
{: .language-python}

4. Client code. The code that manipulate objects in the hierarchy using
the interface defined by Component. 

~~~
def main():

    argon = Sphere()
    krypton = Sphere()
    liquid_crystal = Gay_Berne()
   
    cluster = Cluster()
    cluster.add_particle(argon)
    cluster.add_particle(krypton)
    cluster.add_particle(liquid_crystal)

    argon.translate()
    krypton.translate()
    cluster.translate()

if __name__ == "__main__":
    main()

~~~
{: .language-python}
