---
title: 'Factory Design Pattern'
teaching: 30
exercises: 0
questions:
- 'How can a method or class defer instantiation to subclasses?'
objectives:
- 'Learn the factory design pattern.'
- 'See an example of the factory design pattern relavant to the Computational Molecular Sciences domain.'
keypoints:
- 'A factory allows writing subclasses that change the way an object
is created.'
- 'Factories promote SOLID design principles, enforcing code to be designed towards an interface instead of towards a specific class.'
---

## Definition

The Factory design pattern is a method that uses a superclass to provide a
common interface for the creation of related objects. Subclasses decide which
object type should be created and of hide the implementation details of such
objects.

## Problem

Consider the adapters we previously designed for the two different trajectory analysis tools: MDAnalysis and MDTraj.
Currently we only have two different adapters built, but who knows how many may be constructed in the future. Modifying code to select just between a few toolkits already starts to get tedious if you have to change it a lot. It woul dbe much easier to just take the name of the toolkit you want to use as a parameter and not have to change any code.
We want to create a single interface that remains constant that encapsulate the creation of the adapters and abstracts the specifics away from any user facing code.

## Solution: the factory design pattern

The factory design pattern suggests to replace direct object construction calls
(i.e. instantiating the MDAnalysisAdapter or MDTrajAdapter classes) with calls to a special class called a factory. This class has a method which returns instances of each adapter and hides the specific instantiation of each.
A factory works best when the different constructed classes work from a shared base class, this allows code to use any object produced by the factory in the same way, without even knowing the type of the object they posses.

We restate the two adapter classes here for easier viewing:
~~~
class MDAnalysisAdapter(TrajectoryAdapter):
	def __init__(self, filename):
		self.trajectory = MDAnalysis.Universe(filename)
		print('Selected MDAnalysis.')
	
	def compute_center_of_mass(self):
		mass_by_frame = np.ndarray(shape=(len(self.trajectory.trajectory), 3))
		for ts in self.trajectory.trajectory:
			mass_by_frame[ts.frame] = self.trajectory.atoms.center_of_mass(compound='segments')
		return mass_by_frame
	
	def compute_radius_of_gyration(self):
		rg_by_frame = np.empty(len(self.trajectory.trajectory))
		for ts in self.trajectory.trajectory:
			rg_by_frame[ts.frame] = self.trajectory.atoms.radius_of_gyration()
		return rg_by_frame
		
		
class MDTrajAdapter(TrajectoryAdapter):
	def __init__(self, filename):
		self.trajectory = md.load_pdb(filename)
		print('Selected MDTraj.')
	
	def compute_center_of_mass(self):
		return 10*md.compute_center_of_mass(self.trajectory)
	
	def compute_radius_of_gyration(self):
		return 10*md.compute_rg(self.trajectory)
~~~
{: .language-python}

First we will want to create a class for out factory:
~~~
class TrajectoryAnalysisFactory():
~~~
{: .language-python}
We simply construct a class with a reasonable name for our factory.
We do not need to create an `__init__` method currently, since we do not have and class attributes, the default `__init__` method is sufficient.

## Error handling

We realize that the LJ and Buckingham potentials have different number of
parameters. For our problem, the LJ potential has two input 
parameters, while the Buckingham potential has three. We can improve the
usability of our library by adding the following checks in the potential
constructors to make sure we input the right number of parameters and we use
the desired keywords:


We can write a similar piece of code for the Buckingham potential.

## Final code

The first line of our function creates a dictionary whose values are classes of
potentials. Note that these have not been instantiated yet. After an if
statements that handles errors, we choose the correct class from the dictionary
using the input arguments to our function. The cls variable contains the class
of potential that we need. The line cls(kwargs) instantiates the class with the
provided arguments (i.e. sigma or epsilon for LJ or A, C and rho for
Buckingham).

One advantage of using dictionaries over if statements is that they can be
coupled with a registry
(http://scottlobdell.me/2015/08/using-decorators-python-automatic-registration/)
and make the construction of the dictionary automatic and improve the
extensibility of the potential library. 