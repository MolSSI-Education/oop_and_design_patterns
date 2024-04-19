# Factory Design Pattern

````{admonition} Overview
:class: overview

Questions:
- How can a method or class defer instantiation to subclasses?

Objectives:
- Learn the factory design pattern.
- See an example of the factory design pattern relavant to the Computational Molecular Sciences domain.
````


## Definition

The Factory design pattern is a method that uses a superclass to provide a
common interface for the creation of related objects. Subclasses decide which
object type should be created and of hide the implementation details of such
objects.

## Problem

Consider the adapters we previously designed for the two different trajectory analysis tools: MDAnalysis and MDTraj.
Currently we only have two different adapters built, but who knows how many may be constructed in the future. Modifying code to select just between a few toolkits already starts to get tedious if you have to change it a lot. It would be much easier to just take the name of the toolkit you want to use as a parameter and not have to change any code.
We want to create a single interface that remains constant that encapsulate the creation of the adapters and abstracts the specifics away from any user facing code.

## Solution: the factory design pattern

The factory design pattern suggests to replace direct object construction calls
(i.e. instantiating the MDAnalysisAdapter or MDTrajAdapter classes) with calls to a special class called a factory. This class has a method which returns instances of each adapter and hides the specific instantiation of each.
A factory works best when the different constructed classes work from a shared base class, this allows code to use any object produced by the factory in the same way, without even knowing the type of the object they posses.

We restate the two adapter classes here for easier viewing:
````{tab-set-code} 

```{code-block} python
class MDAnalysisAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = mda.Universe(filename)
        print('Selected MDAnalysis')

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
```
````


````{tab-set-code} 

```{code-block} python
class MDTrajAdapter(TrajectoryAdapter):
    def __init__(self, filename):
        self.trajectory = mdt.load_pdb(filename)
        print('Selected MDTraj.')
	
    def compute_center_of_mass(self):
        return 10 * mdt.compute_center_of_mass(self.trajectory)
	
    def compute_radius_of_gyration(self):
        return 10 * mdt.compute_rg(self.trajectory)
```
````


First we will want to create a new module for our factory, called `factory.py`.
We need to first import our adapter abstract class so our factory knows which type of objects it is building.

````{tab-set-code} 

```{code-block} python
from trajectory_adapter import TrajectoryAdapter
from mdanalysis_adapter import MDAnalysisAdapter
from mdtraj_adapter import MDTrajAdapter
```
````


Then we create a factory class that will act as our factory. This class will hold a few functions going forward. First and foremost is the factory method, which will take in a toolkit name (the library we want to use) and a filename to initialize an adapter.

````{tab-set-code} 

```{code-block} python
class MDFactory:
    def md_factory(self, md_toolkit, filename):
```
````

There are a number of ways we can handle this method, the first way will simply use a series of if statements to determine which adapter needs to be used.
````{tab-set-code} 

```{code-block} python
class MDFactory:
    def md_factory(self, md_toolkit, filename):
        if md_toolkit == 'MDTraj':
            return MDTrajAdapter(filename)
        elif md_toolkit == 'MDAnalysis':
            return MDAnalysisAdapter(filename)
        else:
            raise TypeError('Toolkit not found')
```
````

An alternative option, that looks a little cleaner is to use dictionaries.
````{tab-set-code} 

```{code-block} python
class MDFactory:
    def md_factory(self, md_toolkit, filename):
        toolkits = {'MDTraj': MDTrajAdapter, 'MDAnalysis': MDAnalysisAdapter}
    
        if md_toolkit not in toolkits.keys():
            raise Exception('Toolkit not found.')
        cls = toolkits[md_toolkit]
        return cls(filename)
```
````


Both of these options will work, however, there is a slight problem. They are very rigid and inflexible. We are using a factory so our code can remain open to further extension, further adapters being written. The current factory needs to be modified every time an adapter is written or sufficiently changed. We can solve this issue through the use of a registry.

We start with the dictionary based factory, and simply move the dictionary to be an empty list, as a class variable, and update our references to it, making our module look like this:
````{tab-set-code} 

```{code-block} python
from trajectory_adapter import TrajectoryAdapter

class MDFactory:
    _toolkits = {}

    def md_factory(self, md_toolkit, filename):
        if md_toolkit not in self._toolkits.keys():
            raise Exception('Toolkit not found.')
        cls = self._toolkits[md_toolkit]
        return cls(filename)
```
````


We can now create a new method that acts on this dictionary to register our toolkits.
````{tab-set-code} 

```{code-block} python
    def register(self, toolkit_name, toolkit_class):
        if not issubclass(toolkit_class, TrajectoryAdapter):
            raise TypeError(f'{toolkit_class} is not a TrajectoryAdapter')
        self._toolkits[toolkit_name] = toolkit_class
```
````

Here we first check that the toolkit being registered is a subclass of our `TrajectoryAdapter`, this ensures that it has the proper methods implemented within it. We then simply add it to the dictionary.

The final step is to register each of our toolkits into the factory as needed.. Within our ascript, we simply create a factory and utilize it for registration and getting an adapter.
````{tab-set-code} 

```{code-block} python
from mdtraj_adapter import MDTrajAdapter
from mdanalysis_adapter import MDAnalysisAdapter
from factory import MDFactory

mdf = MDFactory()
mdf.register('MDTraj', MDTrajAdapter)
mdf.register('MDAnalysis', MDAnalysisAdapter)

md = mdf.md_factory('MDTraj', 'protein.pdb')
print(f'Center of mass:\n{md.compute_center_of_mass()}')
print(f'Radius of Gyration:\n{md.compute_radius_of_gyration()}')
```
````


Now any time a new adapter is created, it can be added to the factory and utilized with minimal impact.
To change between which adapter is being used, the only change that needs to be made to the script is the toolkit name:

````{tab-set-code} 

```{code-block} python
#md = mdf.md_factory('MDTraj', 'protein.pdb')
md = mdf.md_factory('MDAnalysis', 'protein.pdb')
```
````
````{admonition} Key Points
:class: key

- A factory allows writing subclasses that change the way an object
````
