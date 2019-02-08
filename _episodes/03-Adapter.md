---
title: "Adapter Design Pattern"
teaching: 30
exercises: 0
questions:
- "How can we utilize pre-existing code without modifying it?"
objectives:
- "Learn the Adapter Design Pattern."
- "See an example of the Adapter Design Pattern relavant to the Computational Molecular Sciences domain."
keypoints:
- "An Adapter allows a class to be used as another interface, without modifying the original class."
- "Adapters promote SOLID design principles, enforcing code to be designed towards an interface instead of towards a specific class."
---
The Adapter design pattern is a method of constructing a wrapper around a pre-existing interface to allow its use by another class.
We can construct an individual adapter for a specific class, or create an abstract interface to base our adapter on.
An abstract interface will allow us to swap out our adapter for any class that inherits from the interface.

An adapter is most useful in two differenct scenarios.
The first is when you would like to re-use pre-existing code that has been rigorously tested without modifying it and thus invalidating the tests.
The second is when you would like to create code that cares more about getting correct results for a set of inputs than how those results are generated.
You may want to create a set of libraries that can generate the results, but you do not want to re-write your main class to suit each new library.
The structure of an Adapter is shown in the figure below.
![adapter](../fig/adapter.png)

As an example, consider the situation where we would like to read a Molecular Dynamics (MD) trajectory from a PDB file and perform some analysis on it.
We currently have an MD library in mind to use, but we may want to use one or more different libraries in the future.
If we do not want to re-write our code for each individual library, we need to come up with a different solution, in this case it will be to create an adapter.
We will use two MD libraries, [MDAnalysis] and [MDTraj] in this example.
We will also be using [NumPy] to structure some of our data.

First, we will import the libraries we want to use.
~~~
import numpy as np
import MDAnalysis
import mdtraj as md
~~~
{: .language-python}

Based on SOLID design principles, we want to code towards an interface, not a specific class, so we want to create an interface to base our adapters on.
We will use Python's abc module to help build the interface, so we will need to import it and build our interface.
~~~
from abc import abstractmethod, ABC

class TrajectoryAdapter(ABC):
	
	@abstractmethod
	def compute_center_of_mass():
		pass
		
	@abstractmethod
	def compute_radius_of_gyration():
		pass
~~~
{: .language-python}
Inheriting from ABC and decorating the methods with `@abstractmethod` ensures that any subclass of TrajectoryAdapter must override both methods.
Any code developed using the listed abstract methods from the interface will now work with any adapter we construct that inherits from `TrajectoryAdapter`.

We will start by building an Adapter that utilizes MDAnalysis.
~~~
class MDAnalysisAdapter(TrajectoryAdapter):
	def __init__(self, filename):
		self.trajectory = MDAnalysis.Universe(filename)
		print('Selected MDAnalysis.')
~~~
{: .language-python}
Here we define the class with a constructor that takes in the name of our PDB file.
We construct a trajectory object by using the `Universe(filename)` method from the MDAnalysis library.
To help determine that we are using the correct adapter, we include a simple print statement informing the User that the MDAnalysis library is selected.

Since both methods in `TrajectoryAdapter` are abstract methods, we must override them.
First we will implement the `compute_center_of_mass` function.
~~~
def compute_center_of_mass(self):
		mass_by_frame = np.ndarray(shape=(len(self.trajectory.trajectory), 3))
		for ts in self.trajectory.trajectory:
			mass_by_frame[ts.frame] = self.trajectory.atoms.center_of_mass(compound='segments')
		return mass_by_frame
~~~
{: .language-python}
We would like to see the result returned as a 2-dimensional array of tuples containing the coordinates of the center of mass in each frame, so we construct our array using NumPy.
We then iterate through the timesteps of our trajectory and calculate the center of mass of the atoms during that frame, adding the result to our array.

We will similarly implement the `compute_radius_of_gyration` function.
~~~
def compute_radius_of_gyration(self):
		rg_by_frame = np.empty(len(self.trajectory.trajectory))
		for ts in self.trajectory.trajectory:
			rg_by_frame[ts.frame] = self.trajectory.atoms.radius_of_gyration()
		return rg_by_frame
~~~
{: .language-python}

We can now use our MDAnalysisAdapter to calculate the center of mass and radius of gyration of a given trajectory.
Using the PDB file provided in the setup section, we can test our code.
~~~
mda = MDAnalysisAdapter('protein.pdb')
print('Center of mass:\n', mda.compute_center_of_mass())
print('Radius of Gyration:\n', mda.compute_radius_of_gyration())
~~~
{: .language-python}
~~~
Selected MDAnalysis.
Center of mass:
 [[60.24883063 51.62894009 28.34133281]
 [60.26522113 51.11509037 27.68827743]
 [60.52356081 50.5282606  27.96596918]
 [60.80518513 49.88031661 26.9061283 ]
 [59.71256198 50.43671635 25.79758075]
 [58.25488384 52.97249614 26.26805148]
 [57.5791976  52.26752831 26.35253427]
 [57.76277669 52.2229901  24.79696899]
 [56.62274735 52.48640321 26.98199884]
 [56.7851466  52.96245182 27.8464323 ]]
Radius of Gyration:
 [24.37682533 23.73724108 23.45717622 23.85568415 23.38966031 21.46212035
 21.95454317 21.96694184 21.22233935 20.48089334]
~~~
{: .output}

We would like to change the library we are using to MDTraj without adjusting the function calls.
First, we will construct another Adapter for MDTraj, implenting both abstract methods.
~~~
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
It is important to notice that the constructor and the method calls have the same definition as those in the MDAnalysisAdapter and perform the same operations, but do so in a different way.
To an external class, the two Adapters are now interchangable.
We will run the same test with the only difference being the adapter we use.
~~~
mda = MDTrajAdapter('protein.pdb')
print('Center of mass:\n', mda.compute_center_of_mass())
print('Radius of Gyration:\n', mda.compute_radius_of_gyration())
~~~
{: .language-python}
~~~
Selected MDTraj.
Center of mass:
 [[60.24882436 51.62893846 28.34134574]
 [60.26521674 51.11507794 27.68828924]
 [60.5235568  50.5282497  27.96598254]
 [60.80518015 49.88031491 26.90614047]
 [59.71255677 50.43671612 25.79759428]
 [58.25487874 52.97249014 26.26806535]
 [57.57919235 52.26752287 26.35254909]
 [57.7627734  52.22298068 24.79698321]
 [56.62274205 52.48639788 26.98201221]
 [56.78514205 52.96243944 27.84644585]]
Radius of Gyration:
 [24.29005389 23.54864649 23.33444906 23.84356938 23.40307752 21.52953619
 21.96178093 21.93913835 21.22483198 20.43491991]
~~~
{: .output}

With adapters for each library, our code is not concerned with how the data it needs is generated, simply that it follows the contract put in place by the interface.

{% include links.md %}

[MDAnalysis]: https://www.mdanalysis.org/
[MDTraj]: http://mdtraj.org/1.9.0/
[NumPy]: http://www.numpy.org/