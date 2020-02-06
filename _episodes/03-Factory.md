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

First we will want to create a new module for our factory, called `factory.py`.
We need to first import our adapter abstract class so our factory knows which type of objects it is building.

~~~
from .trajectoryadapter import TrajectoryAdapter
~~~
{: .language-python}

Then we create a funciton that will act as our factory. We want this function to take in arguments for our adapters, along with a label so it knows which adapter to create.

~~~
def trajectory_factory(trajectory_toolkit, **kwargs):
~~~
{: .language-python}
There are a number of ways we can handle this method, the first way will simply use a series of if statements to determine which adapter needs to be used.
~~~
def trajectory_factory(trajectory_toolkit, **kwargs):
    if trajectory_toolkit == 'MDTraj':
        traj_analysis = MDTrajAdapter(kwargs['filename'])
    elif trajectory_toolkit == 'MDAnalysis':
        traj_analysis = MDAnalysisAdapter(kwargs['filename'])
    else:
        raise TypeError('Toolkit not found')
~~~
{: .language-python}
An alternative option, that looks a little cleaner is to use dictionaries.
~~~
def trajectory_factory(trajectory_toolkit, **kwargs):
    traj_toolkits = {'MDTraj': MDTrajAdapter, 'MDAnalysis': MDAnalysisAdapter}
    
    if trajectory_toolkit not in traj_toolkits.keys():
        raise TypeError('Toolkit not found')
    
    traj_analysis = traj_toolkits[trajectory_toolkit](kwargs['filename'])
    
    return traj_analysis
~~~
{: .language-python}

Both of these options will work, however, there is a slight problem. They are very rigid and inflexible. We are using a factory so our code can remain open to further extension, further adapters being written. The current factory needs to be modified every time an adapter is written or sufficiently changed. We can solve this issue through the use of a registry.

We start with the dictionary based factory, and simply move the dictionary to be an empty list, outside of the function, and update our references to it, making our module look like this:
~~~
from .trajectoryadapter import TrajectoryAdapter

_toolkits = {}

def trajectory_factory(trajectory_toolkit, **kwargs):
    traj_toolkits = {'MDTraj': MDTrajAdapter, 'MDAnalysis': MDAnalysisAdapter}
    
    if trajectory_toolkit not in _toolkits.keys():
       raise TypeError('Toolkit not found')
    
    traj_analysis = _toolkits[trajectory_toolkit](kwargs['filename'])
    
    return traj_analysis
~~~
{: .language-python}

We can now create a new method that acts on this dictionary to register our toolkits.
~~~
def register(toolkit_name, toolkit_class):
    if not issubclass(toolkit_class, TrajectoryAdapter):
        raise TypeError('{0} is not a TrajectoryAdapter'.format(toolkit_class))
    _toolkits[toolkit_name] = toolkit_class
~~~
{: .language-python}
Here we first check that the toolkit being registered is a subclass of our `TrajectoryAdapter`, this ensures that it has the proper methods implemented within it. We then simply add it to the dictionary.

The final step is to make each of our current toolkits register themselves into the factory. Within our adapter module, we simply import the `register` method and use it for each adapter.
~~~
from .factory import register

register('MDTraj', MDTrajAdapter)
register('MDAnalysis', MDAnalysisAdapter)
~~~
{: .language-python}

Now any time a new adapter is created, it can be added to the factory without any modifications outside of its own module, it just imports the register function and utilizes it.

## Final code
The final code for the `factory` module is included here for clarity.
~~~
from .trajectoryadapter import TrajectoryAdapter

_toolkits = {}

def register(toolkit_name, toolkit_class):
    if not issubclass(toolkit_class, TrajectoryAdapter):
        raise TypeError('{0} is not a TrajectoryAdapter'.format(toolkit_class))
    _toolkits[toolkit_name] = toolkit_class

def trajectory_factory(trajectory_toolkit, **kwargs):
    if trajectory_toolkit not in _toolkits.keys():
       raise TypeError('Toolkit not found')
       
    traj_analysis = _toolkits[trajectory_toolkit](kwargs['filename'])
    
    return traj_analysis
~~~
{: .language-python}
