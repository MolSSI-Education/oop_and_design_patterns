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

Imagine you are creating a Monte Carlo or molecular dynamics library. An
important component of such library is the pairwise energy computation. You
might be interested in implementing at least two potentials: the 12-6
Lennard-Jones potential and the Buckingham potential.

You implement such functional forms using the following classes:


~~~
class LJ:
   def __init__(self, epsilon, sigma):
       self.sigma = sigma
       self.epsilon = epsilon

   def get_energy(self, r):
       return 4 * self.epsilon * ((self.sigma / r)**12 - (self.sigma / r)**6)


class Buckingham:
   def __init__(self, rho, A, C):
       self.rho = rho
       self.A = A
       self.C = C

   def get_energy(self, r):
       return self.A * np.exp(-r / self.rho) - self.C / r**6

~~~
{: .language-python}

You quickly realize that you should design your code such that is extensible to
many different functional forms. Although you could have a class that defines
each potential and have the client (i.e. code that uses your library)
instantiate each class every time, that might be cumbersome in the long term.
It would be nice to have a common interface that encapsulates the creation of
various potential functional forms so that your code has a consistent interface
over time.

## Solution: the factory design pattern

The factory design pattern suggests to replace direct object construction calls
(i.e. instantiating the LJ or Buckingham classes) with calls to a special
method called a factory. This method returns instances of each potential class
and hides the implementation of the potential classes.
	
The factory method for our simulation potentials could look like:

~~~
def potential_factory(potential_type, **kwargs):
    if potential_type == 'LJ':
        return LJ(**kwargs)

    elif potential_type == 'Buckingham':
        return Buckingham(**kwargs)

    else:
        raise Exception('Potential type not found')
~~~
{: .language-python}

In this way, our client code can call our factory class as follows:

~~~
buckingham_potential = potential_factory('Buckingham', A=4.0, rho=10.0, C=10)
energy = buckingham_potential.get_energy(r=10.0)
~~~
{: .language-python}

We have created a common interface (the factory  method) that returns potential
objects given a chosen keyword (i.e. Buckingham or LJ).

We realize that all potential functional classes must have a `get_energy`
method that computes the energy for a given distance. Without this method, our factory would
return inconsistent objects. We can achieve this by using an abstract base
class that forces the implementation of an abstract method to its children
classes


~~~
from abc import ABC, abstractmethod


class Potential(ABC):
   @abstractmethod
   def get_energy(self):
       pass
~~~
{: .language-python}

If we make our LJ and Buckingham  classes children of the Potential class, we
force them to implement the get_energy method and thus our factory would return
consistent objects.  We are done with implementing our factory!

## Error handling

We realize that the LJ and Buckingham potentials have different number of
parameters. For our problem, the LJ potential has two input 
parameters, while the Buckingham potential has three. We can improve the
usability of our library by adding the following checks in the potential
constructors to make sure we input the right number of parameters and we use
the desired keywords:

~~~
for key in kwargs:
    if key not in ['epsilon', 'sigma']:
        raise KeyError('LJ potential: Must input epsilon and sigma')
~~~
{: .language-python}

We can write a similar piece of code for the Buckingham potential.

## Final code

The final version of our code is


~~~

from abc import ABC, abstractmethod
import numpy as np


class Potential(ABC):
   @abstractmethod
   def get_energy(self):
       pass


class LJ(Potential):
   def __init__(self, **kwargs):
       for key in kwargs:
           if key not in ['epsilon', 'sigma']:
               raise KeyError('LJ potential: Must input epsilon and sigma')

       self.sigma = kwargs['sigma']
       self.epsilon = kwargs['epsilon']

   def get_energy(self, r):
       return 4 * self.epsilon * ((self.sigma / r)**12 - (self.sigma / r)**6)


class Buckingham(Potential):
   def __init__(self, **kwargs):
       for key in kwargs:
           if key not in ['A', 'C', 'rho']:
               raise KeyError('Buckingham potential: Must input A, C and rho')

       self.rho = kwargs['rho']
       self.A = kwargs['A']
       self.C = kwargs['C']

   def get_energy(self, r):
       return self.A * np.exp(-r / self.rho) - self.C / r**6


def potential_factory(potential_type, **kwargs):
    if potential_type == 'LJ':
        return LJ(**kwargs)

    elif potential_type == 'Buckingham':
        return Buckingham(**kwargs)

    else:
        raise Exception('Potential type not found')

# Client code example below

buck_potential = potential_factory('Buckingham', A=4.0, rho=10.0, C=10)
energy = buck_potential.get_energy(r=10.0)

~~~
{: .language-python}

## Alternative implementation of the factory method

Instead of using a set of if statements, we could use a Python dictionary to implement the same functionality as follows

~~~
def potential_factory(potential_type, **kwargs):
   cls_dict = dict(LJ=LJ, Buckingham=Buckingham)

   if potential_type not in cls_dict.keys():
       raise Exception('Potential type not found')

   cls = cls_dict[potential_type]

   cls_instance = cls(**kwargs)
   return cls_instance
~~~
{: .language-python}

The first line of our function creates a dictionary whose values are classes of
potentials. Note that these have not been instantiated yet. After an if
statements that handles errors, we choose the correct class from the dictionary
using the input arguments to our function. The cls variable contains the class
of potential that we need. The line cls(**kwargs) instantiates the class with the
provided arguments (i.e. sigma or epsilon for LJ or A, C and rho for
Buckingham).

One advantage of using dictionaries over if statements is that they can be
[coupled with a registry]
(http://scottlobdell.me/2015/08/using-decorators-python-automatic-registration/)
and make the construction of the dictionary automatic and improve the
extensibility of the potential library. 
