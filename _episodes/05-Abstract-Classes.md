# Abstract Classes

````{admonition} Overview
:class: overview

Questions:
- What are Abstract Classes?
- Why should I use Abstract Classes?

Objectives:
- Understand the concepts behind Abstract Classes.
````

## Abstract Classes
An abstract class is a way to define how a class will be designed while allowing some or all of its methods to remain unimplemented.
An abstract class can be inherited by another abstract class, in which case it is usually extended into a larger abstract class by the child abstract class, or it can be inherited by a class, which implements at least all of the abstract methods from the abstract class.
An abstract class creates a structure for code to look for without defining how the methods are implemented.
Any methods implemented within an abstract class can be used by their children or overridden to achieve new behavior.

In python, we construct abstract classes using a package from the standard library.
We first need to import the abstract base class and abstract method decorator.
````{tab-set-code} 

```{code-block} python
from abc import ABC, abstractmethod
```
````

`ABC` is the abstract base class.
Inheriting from the abstract base class enforces that any child classes of the abstract class must have some implementation of any abstract methods.
````{tab-set-code} 

```{code-block} python
class AbstractClass(ABC):
```
````

Here we define the name of our abstract class.
Note `ABC` in the class definition to denote the inheritance.
Then we create a definition for each method we would like any of the abstract classes children to implement and any methods we would like to pass on.
````{tab-set-code} 

```{code-block} python
    @abstractmethod
    def first_method(self):
        pass

    @abstractmethod
    def second_method(self):
        pass
        
    def third_method(self):
        print('This method is implemented for use by any child classes.')
```
````

Since we have used the `@abstractmethod` decorator for some of these methods, any child class must implement a method with the same name as the abstract class.

If we try and create an instance of this object, we will notice a problem.
````{tab-set-code} 

```{code-block} python
abstract_class = AbstractClass()
```
````

````{tab-set-code} 

```{code-block} python
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: Can't instantiate abstract class AbstractClass with abstract methods first_method, second_method
```
````
Since there are unimplemented methods, we canot create the instance of the object.

A child class that inherits our abstract class will need to override `first_method` and `second_method` if it is not itself an abstract class.

````{tab-set-code} 

```{code-block} python
class ChildClass(AbstractClass):
    def first_method(self):
        print('Implementing the first method inherited from my parent.')
    
    def second_method(self):
        print('Implementing the second method inherited from my parent.')
```
````


Note that we did not implement the `third_method` in the child class. Since `third_method` is not an abstract method, we can utilize the implementation provided by our parent.

Lets try and create and instance of the child class.
````{tab-set-code} 

```{code-block} python
child_class = ChildClass()
```
````


We should see a blank output as the object has been created correctly.
We can now use all three methods on our child object.
````{tab-set-code} 

```{code-block} python
child_class.first_method()
child_class.second_method()
child_class.third_method()
```
````

````{tab-set-code} 

```{code-block} output
Implementing the first method inherited from my parent.
Implementing the second method inherited from my parent.
This method is implemented for use by any child classes.
```
````


If we wish to override the method provided by the abstract class, we can simply re-write the method.
````{tab-set-code} 

```{code-block} python
class ChildClass(AbstractClass):
    def first_method(self):
        print('Implementing the first method inherited from my parent.')
    
    def second_method(self):
        print('Implementing the second method inherited from my parent.')
    
    def third_method(self):
        print('Overriding the third method provided by my parent.')
```
````


Re-instantiating the object and running the three methods, we will see the overridden output.
````{tab-set-code} 

```{code-block} python
child_class = ChildClass()

child_class.first_method()
child_class.second_method()
child_class.third_method()
```
````

````{tab-set-code} 

```{code-block} output
Implementing the first method inherited from my parent.
Implementing the second method inherited from my parent.
Overriding the third method provided by my parent.
```
````


## Molecule Example

Given that we have already created a `Molecule` class, let us assume that we wanted it to inherit from a base class.
Consider the following base class as an example:
````{tab-set-code} 

```{code-block} python
class MolBase(ABC):
    @abstractmethod
    def compute_center_of_mass():
        pass
```
````


`MolBase` inherits from the abstract base class and has a single abstract method.
Trying to instantiate an object of type `MolBase`, we get an error, as `MolBase` has an abstract method.
````{tab-set-code} 

```{code-block} python
mol1 = MolBase()
```
````

````{tab-set-code} 

```{code-block} output
TypeError                                 Traceback (most recent call last)
/tmp/ipykernel_1250/3850661068.py in <module>
----> 1 mol1 = MolBase()

TypeError: Can't instantiate abstract class MolBase with abstract methods compute_center_of_mass
```
````


If we change the `Molecule` class to inherit from `MolBase` and try and create a molecule we will get a similar error:
````{tab-set-code} 

```{code-block} python
class Molecule(MolBase):
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates

    @property
    def symbols(self):
        return self._symbols
        
    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols
        self._update_num_atoms()

    def _update_num_atoms(self):
        self.num_atoms = len(self.symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}\nnumber of atoms: {self.num_atoms}'

mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
print(mol1)
```
````

````{tab-set-code} 

```{code-block} output
TypeError                                 Traceback (most recent call last)
/tmp/ipykernel_1250/1044270064.py in <module>
----> 1 mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
      2 print(mol1)

TypeError: Can't instantiate abstract class Molecule with abstract methods compute_center_of_mass
```
````

To implement the center of mass for `Molecule`, we will also need to pass in a set of masses for the atoms, which we will do in the form of a list.

The full class with masses and an implemented `compute_center_of_mass()` method:
````{tab-set-code} 

```{code-block} python
class Molecule(MolBase):
    def __init__(self, name, charge, symbols, coordinates, masses):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
        self.masses = masses

    @property
    def symbols(self):
        return self._symbols
        
    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols
        self._update_num_atoms()

    def _update_num_atoms(self):
        self.num_atoms = len(self.symbols)

    def compute_center_of_mass(self):
        com_x = 0
        com_y = 0
        com_z = 0
        total_mass = 0
        for mass, coordinate in zip(self.masses, self.coordinates):
            total_mass += mass
            com_x += mass*coordinate[0]
            com_y += mass*coordinate[1]
            com_z += mass*coordinate[2]
        return [com_x/total_mass, com_y/total_mass, com_z/total_mass]

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}\nnumber of atoms: {self.num_atoms}'
```
````

We can create an instance of `Molecule` and try and run the method.
````{tab-set-code} 

```{code-block} python
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]], masses=[15.999, 1.00784, 1.00784])
print(mol1)
print(mol1.compute_center_of_mass())
```
````

````{tab-set-code} 

```{code-block} output
name: water molecule
charge: 0.0
symbols: ['O', 'H', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
number of atoms: 3
[0.0, 0.05594548446045114, 0.05594548446045114]
```
````


Having implemented the abstract method, our class is instantiable.
Note however, that python does not concern itself with the specifics of the implementation.
Python's only concern is that a non-abstract method exists with the same name.
The input parameters and implementation will not be checked.
````{admonition} Key Points
:class: key

- Abstract Classes
- Abstract Methods
- Overwritting Abstract Methods
````
