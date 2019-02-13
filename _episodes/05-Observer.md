---
title: 'Observer Design Pattern'
teaching: 30
exercises: 0
questions:
- 'How can an object notify its own state to an open-ended number of objects?'
objectives:
- 'Learn the observer design pattern.'
- 'See an example of the observer design pattern relavant to the Computational Molecular Sciences domain.'
keypoints:
- 'The observer design pattern provides a way for the subject to notify an open-ended number of objects about its own state'
---

## Definition

Define a one-to-many dependency between objects so that when one object changes
state, all its dependents are notified and updated automatically.

## Problem

Imagine you are creating a classical Monte Carlo library. You might
be interested in retrieving information such as Cartesian coordinates
snapshots, potential and kinetic energies, pressures, densities and
other thermodynamic quantities at any given moment during the simulation. 
Aditionally, you would like to have 
a checkpointing mechanism that captures the current states of the simulation
in case it crashes.

In general, you would like to have many ways of looking at the current
state of your system. Your code should make it easy to incorporate new
indicators of the state of the system.

## Solution: the observer design pattern

The observer design pattern provides a way for your simulation system (the subject) to notify its attached observers about its current state. All
observers are notified whenever the subject undergoes a change in state.

The subject knows about its observers and provides an interface for
attaching and detaching observers. We can create an abstract base clase
for defining all subjects as:

~~~
class Subject:
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        if observer not in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        if observer in self._observers:
            self._observers.remove(observer)

    def notify(self):
        for observer in self._observers:
            observer.update(self)
~~~
{: .language-python}

Now, we need the actual concrete subject that our observers will look at. For our
case, it could be our simulation system (i.e. all particles). Thus, we need to
define a concrete subject that inherits from our abstract subject class

~~~
class System(Subject):

    def __init__(self, mol_number=0):
        Subject.__init__(self)

	# Initialize particle positions
        self._xyz = np.zeros((mol_number, 3))

        self.current_step = 0

	# Compute initial energy
        self.compute_energy()

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz
        self.current_step += 1
        self.compute_energy()
        self.notify()

    def compute_energy(self):
        # Ideally we would implement a full energy computation
        self.energy = np.random.rand()
~~~
{: .language-python}


Additionally, we need a base class for all observers. It should require a
function that retrieves what the observer knows about the subject.

~~~
class Observer(ABC):

    @abstractmethod
    def __init__(self, print_freq=1):
        pass

    @abstractmethod
    def update(self, subject):
        pass
~~~
{: .language-python}

And now we can implement two concrete observers: one that retrieves the Cartesian coordinates
and another that gets the energy.

~~~
class XYZ_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing xyz subject.xyz')


class energy_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing energy', subject.energy)

~~~
{: .language-python}


Now, you could use the above classes in your Monte Carlo client code
like this

~~~
# Create the system with a desired number of particles
mc_system = System(mol_number=15)

xyz_observer = XYZ_observer(print_freq=2)
energy_observer = energy_observer()

mc_system.attach(xyz_observer)
mc_system.attach(energy_observer)

total_mcsteps = 5
for this_mcstep in range(total_mcsteps):

    # Ideally you would perform some MC moves here
    # that alter the state of the system and
    # don't violate detailed balance!
    mc_system.xyz += 0.1
~~~
{: .language-python}

## Final code

The final version of our code is

~~~

from abc import ABC, abstractmethod
import numpy as np


class Subject:
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        if observer not in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        if observer in self._observers:
            self._observers.remove(observer)

    def notify(self):
        for observer in self._observers:
            observer.update(self)


class Observer(ABC):

    @abstractmethod
    def __init__(self, print_freq=1):
        pass

    @abstractmethod
    def update(self, subject):
        pass


class System(Subject):

    def __init__(self, mol_number=0):
        Subject.__init__(self)
        self._xyz = np.zeros((mol_number, 3))
        self.current_step = 0
        self.compute_energy()

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz
        self.current_step += 1
        self.compute_energy()
        self.notify()

    def compute_energy(self):
        # Ideally we would implement a full energy computation
        self.energy = np.random.rand()


class XYZ_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing xyz subject.xyz')


class energy_observer(Observer):

    def __init__(self, print_freq=1):
        self.print_freq = print_freq

    def update(self, subject):

        if subject.current_step % self.print_freq == 0:
            print('Printing energy', subject.energy)


# Create the system with a desired number of particles
mc_system = System(mol_number=15)

xyz_observer = XYZ_observer(print_freq=2)
energy_observer = energy_observer()

mc_system.attach(xyz_observer)
mc_system.attach(energy_observer)

total_mcsteps = 5
for this_mcstep in range(total_mcsteps):

    # Ideally you would perform some MC moves here
    # that alter the state of the system and
    # don't violate detailed balance!
    mc_system.xyz += 0.1


~~~
{: .language-python}

