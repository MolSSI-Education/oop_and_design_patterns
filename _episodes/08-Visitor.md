---
title: 'Visitor Design Pattern'
teaching: 30
exercises: 0
questions:
- 'How to separate the methods from the objects that use them?'
keypoints:
- 'The visitor design pattern lets you define operations without changing the classes of the elements that use such operations.'
objectives:
- 'Understand the visitor design pattern.'
- 'Show an example of the visitor design pattern relavant to the Computational Molecular Sciences domain.' 
---

## Problem

This example will build off of the composite [design pattern example] (https://molssi-education.github.io/oop_and_design_patterns/07-Composite/index.html). 

Imagine you would like to add a solid body rotation Monte Carlo move to your program.
As you start modifying your code, you'll realize that you need to modify your Component
base class and each of the Patchy, Cluster and other children classes:

~~~
from abc import ABC, abstractmethod

class Component(ABC):

    '''Base class for a Component. This is the 'component'
    base class in the pattern.'''

    @abstractmethod
    def translate(self):
        pass

    @abstractmethod
    def rotate(self):
	pass

class Patchy(Component):

    '''A patchy particle is a type of particle.'''

    def translate(self):
        print('translating sphere')

    def rotate(self):
	pass

class Cluster(Component):

    '''A patchy particle is a type of particle.'''

    def __init__(self):
        self._particles = []

    def translate(self):
        for children in self._particles:
            children.translate()

    def rotate(self)
	pass

    def add_particle(self, particle):
        self._particles.append(particle)

~~~
{: .language-python}

Each time we need to add a new move, we will need to change all the existing Particle
classes to add that new move. This might become unpractical the more derived classes you have.

## Solution: the visitor design pattern

The visitor pattern lets you define a new operation without changing the
classes of the elements on which it operates. To implent it, we need
the following elements

1. Visitor: this is an abstract class that is used to declare the methods 
that will be used by the visitable clases. In our example, the Visitor abstract class
will be an abstract Move_Type class that will contain the abstract methods for translation
and rotation. 

~~~
class Move_Type(ABC):

    @abstractmethod
    def sphere(self):
        pass

    @abstractmethod
    def gay_berne(self):
        pass

    @abstractmethod
    def cluster(self):
        pass
~~~
{: .language-python}


2. Concrete visitor: This implements the abstract class. We need to define a
concrete visitor for each type of visitor. For instance, a concrete visitor would be
a "Translator" class, which will define the translation methods for each type
of particle, as defined by its parent class.

~~~
class Translator(Move_Type):
    def sphere(self, sphere):
        print('Translating sphere should be very straightforward')

    def gay_berne(self, gay_berne):
        print('Translating solid should involve COM computation and subsequent translation')

    def cluster(self, cluster):
        print('Translating cluster should involve COM computation and subsequent translation')
        for particle in cluster.particles:
            print('Translating particle', particle) 
~~~
{: .language-python}

For a rotation move, we need to create another class that implements the rotation
routines for each type of particle. 

~~~
class Rotator(Move_Type):
    def sphere(self, sphere):
        print('Rotating sphere is nonsense due to rotation invariance')

    def gay_berne(self, gay_berne):
        print('Rotating Gay Berne using Eulerian angles', gay_berne)

    def cluster(self, cluster):
        print('Rotating cluster using Eulerian angles')
        for particle in cluster.particles:
            print('Rotating particle', particle) 
~~~
{: .language-python}

Adding a new move would simply involve creating a new class and, for each type of particle, 
defining the method that would correctly handle the transformation for that particular
type. If a particular move does not make sense for a type of particle, one could
simply not implement it, raise an exception, or do whatever you want (for instance
the rotation move for a sphere).

3. Visitable: is an abstraction which defines the method that will allow a children
object to be visited. In our example this would simply be an abtract Particle
class that defines a generic transform method.

~~~
class Particle(ABC):

    '''Base class for a Particle. Particles
    can be translated, rotated, reflected, etc'''

    @abstractmethod
    def transform(self, transformer):
        pass
~~~
{: .language-python}

4. Concrete Visitable: A class that implement the method defined by Visitable. Note
that the visitor object is an argument that is passed to the generic function transform. This 
generic function has to be implemented in each particle type. 

The function move requires as argument an instance of a class Particle, as required
by the abstract Move_Type class. For instance, 
instance, the Translator.sphere method needs an instance
of the Sphere class. This instance is provided in Sphere.move class by the object self, which
represent the current object instance.

~~~
class Sphere(Particle):

    '''A Sphere is a type of particle, for instance, 
    point particles with Lennard-Jones potential'''

    def move(self, transformer):
        transformer.sphere(self) 

class Gay_Berne(Particle):

    '''A Gay Berne particle is an anisotropic model
    of the 12-6 Lennard-Jones potential. It is used
    exstensively to model liquid crystals'''

    def move(self, transformer):
        transformer.gay_berne(self) 
~~~
{: .language-python}

5. Object Structure. This is a class containing all the objects that can be
   visited. It offers a mechanism to iterate through all the elements. 
In can be a complex structure, such
as a composite object

~~~
class Cluster(Particle):

    '''A cluster is a collection of particles. It can
    be composed of spherical, anisotropic, patchy or 
    any other particle'''

    def __init__(self):
        self.particles = []

    def move(self, transformer):
        transformer.cluster(self) 

    def add_particle(self, particle):
        self.particles.append(particle)
~~~
{: .language-python}

Finally, the client code that uses the aforementioned classes is

~~~
def main():

    translate = Translator()
    rotate = Rotator()

    argon = Sphere()
    krypton = Sphere()
    liquid_crystal = Gay_Berne()

    cluster = Cluster()
    cluster.add_particle(argon)
    cluster.add_particle(krypton)

    argon.transform(translate)
    argon.transform(rotate)

    krypton.transform(translate)
    krypton.transform(rotate)

    liquid_crystal.transform(translate)
    liquid_crystal.transform(rotate)

    cluster.transform(translate)
    cluster.transform(rotate)

if __name__ == "__main__":
    main()
~~~
{: .language-python}

## Implementation in Python

As a dynamic language, Python offers functionality that ends up simplifying
and sometimes nullyfing the need for design patterns as described in the
classic book Design Patterns: Elements of Object-Oriented Software. For instance, the visitor pattern can be implemented in Python as

~~~
from abc import ABC, abstractmethod

class Particle:

    '''Base class for a Particle. Particles
    can be translated, rotated, reflected, etc'''

    def transform(self, transformer):
        method_name = 'transform_{}'.format(self.__class__.__name__.lower())
        try:
            visit = getattr(transformer, method_name)
        except AttributeError:

            print("WARNING: The transformer {} does not contain the method {} for type {}".format(transformer.__class__.__name__, method_name, self.__class__.__name__))
            return
        return visit(self)

class Sphere(Particle):

    '''A Sphere is a type of particle, for instance, 
    point particles with Lennard-Jones potential'''
    def __init__(self, center):
        self.center = center

class Gay_Berne(Particle):

    '''A Gay Berne particle is an anisotropic model
    of the 12-6 Lennard-Jones potential. It is used
    exstensively to model liquid crystals'''
    def __init__(self, center, vector):
        self.center = center 
        self.vector = vector 

class Cluster(Particle):

    '''A cluster is a collection of particles. It can
    be composed of spherical, anisotropic, patchy or 
    any other particle'''

    def __init__(self):
        self.particles = []

    def add_particle(self, particle):
        self.particles.append(particle)

class Translator:
    def transform_sphere(self, sphere):
        print('Translating sphere should be very straightforward')
        sphere.center[0] += 1.0

    def transform_gay_berne(self, gay_berne):
        print('Translating solid should involve COM computation and subsequent translation')

    def transform_cluster(self, cluster):
        print('Translating cluster should involve COM computation and subsequent translation')
        for particle in cluster.particles:
            particle.transform(self) 


class Rotator:

    def transform_gay_berne(self, gay_berne):
        print('Rotating Gay Berne using Eulerian angles')
        gay_berne.vector[0] = 1

    def transform_cluster(self, cluster):
        print('Rotating cluster using Eulerian angles')
        for particle in cluster.particles:
            particle.transform(self)

def main():

    translator = Translator()
    rotator = Rotator()

    origin = [0.0, 0.0]
    vector = [1.0, 0.0]

    argon = Sphere(origin)
    krypton = Sphere(origin)
    liquid_crystal = Gay_Berne(origin, vector)
   
    cluster = Cluster()
    cluster.add_particle(argon)
    cluster.add_particle(krypton)
    cluster.add_particle(liquid_crystal)

    argon.transform(translator)
    argon.transform(rotator)
    
    krypton.transform(translator)
    krypton.transform(rotator)

    liquid_crystal.transform(translator)
    liquid_crystal.transform(rotator)

    cluster.transform(translator)
    cluster.transform(rotator)

if __name__ == "__main__":
    main()

~~~
{: .language-python}
