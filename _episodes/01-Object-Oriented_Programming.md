---
title: "Object Oriented Programming"
teaching: 0
exercises: 0
questions:
- "What is Object Oriented Programming?"
- "Why should I use Object Oriented Programming?"
objectives:
- "Understand the concepts behind Object Oriented Programming."
keypoints:
- "Encapsulation"
- "Inheritance"
---

{% include links.md %}
## Object Oriented Programming
Object Oriented Programming (OOP) is a style of programming that promotes the creation of objects that contain data and methods that act on their data.
A program developed using OOP will consist of a number of objects that interact with one another.
The primary concepts of OOP are Encapsulation, Data Abstraction, Polymorphism, and Inheritance.
This introduction will briefly cover each of these concepts and why they are useful to a software developer.


## Encapsulation
Encapsulation is the concept of enclosing the data and methods within a class.
A class will consist of a set of data (variables) and a set of methods that interact with the data.
There are a number of benefits to creating classes.
1. It aids in the understanding of the code being developed.
It is much easier to understand how an object will behave if all the data and methods that interact with that data are enclosed within its class.
It means all of the information about that object is grouped into a single location.
2. As an extension of the first benefit, anyone developing code that utilizes your objects will have a clear understanding of how they are allowed to interact with it.
The methods they can use will be located within the class of the object.
3. It promotes security of data by restricting the ways to access data to the specific methods within the class.

Classes are used in code to provide a general definition of an object.
We call an object an instance of a class, meaning it takes the structure of the class and fills in any necessary data.
Using classes as a general definition allows you to construct many instances of the class with little work.
Consider the following code to define a molecule without using a class.
<!---
TODO: Develop sample of code creating multiple objects without defining classes.
--->

~~~
molecule_name = "water dimer"
molecule_charge = 0.0
molecule_symbols = ["O", "H", "H"]

print('name: ', molecule_name, '\ncharge:', molecule_charge, '\nsymbols:', molecule_symbols)

molecule2_name = "He"
molecule2_charge = 0.0
molecule2_symbols = ["He"]

print('name: ', molecule2_name, '\ncharge:', molecule2_charge, '\nsymbols:', molecule2_symbols)
~~~
{: .language-python}
For each new molecule we want to build using this method, we need to create a new variable name for each of the variables and redefine how we are printing them.
For a single instance, this may not seem terrible, but what if we need to make tens or hundreds of instances of the same type of object?
The amount of additional code that needs to be written grows very quickly.

Now consider code that creates and instantiates a class multiple times instead.
<!---
TODO: Develop sample class to create the multiple objects instead.
--->

~~~
class Molecule:
	def __init__(self, name, charge, symbols):
		self.name = name
		self.charge = charge
		self.symbols = symbols
		
	def print_molecule(self):
		print('name: ', self.name, '\ncharge:', self.charge, '\nsymbols:', self.symbols)
		
mol1 = Molecule('water dimer', 0.0, ["O", "H", "H"])
mol2 = Molecule('He', 0.0, ["He"])

mol1.print_molecule()
mol2.print_molecule()
~~~
{: .language-python}

The construction and use of an object constructed through a class is simpler and more intuitive that trying to construct one from an arbitrary set of variables.
Now anytime we wish to create a new molecule object and print out its values, we only require two new lines of code.

Stepping through the class example, we will briefly cover how to construct a class and instantiate it as an object.
The first line
~~~
class Molecule:
~~~
{: .language-python}
is defining the name of the class as Molecule.
We then have two methods, the first one is called a constructor and it is called whenever you are instantiating an object of the class.
We have the definition of the constructor
~~~
	def __init__(self, name, charge, symbols):
~~~
{: .language-python}
that has three parameters: name, charge, and symbols.
The parameters of a constructor are required anytime you want to create an instance of the class.
These can have default values if they are non-required.
The next three lines
~~~
		self.name = name
		self.charge = charge
		self.symbols = symbols
~~~
{: .language-python}
set the value of the local object variables to the value of the parameters.

The second method is simply printing our Molecule object, similar to how we defined it without any class, but it now will work for each instance of a Molecule without any modification.
~~~
	def print_molecule(self):
		print('name: ', self.name, '\ncharge:', self.charge, '\nsymbols:', self.symbols)
~~~
{: .language-python}

Finally, we need to create the instances of our Molecule class and call their print methods.
~~~
mol1 = Molecule('water dimer', 0.0, ["O", "H", "H"])
mol2 = Molecule('He', 0.0, ["He"])

mol1.print_molecule()
mol2.print_molecule()
~~~
{: .language-python}


## Data Abstraction

<!--- Private vs Public --->
Data abstraction is the concept of hiding variables behind methods to simplify code, increase development efficiency and improve security.


Generally, variables can be public, in which case they are directly modifiable, or private, meaning interaction with their values is only possible through internal class methods.
In python, there are no explicitly public or private variables, all variables are accessible within an object.
However, the predominantly accepted practice is to treate names prefixed with an underscore as non-public. (See [Python Private Variables](https://docs.python.org/3/tutorial/classes.html#private-variables))

<!--- TODO come up with a code example here --->
<!--- TODO Do something to show the differences in interacting with variables through the variable itself and through getters/setters. --->

## Inheritance
Inheritance is the principle of extending a class to add capabilities without modifying the original class.
We call the class that is being inherited the parent, and the class that is inheriting the child.
The child class obtains the properties and behaviors of its parent unless it overrides them.

In coding terms, this means a class that inherits from a parent class by default will contain all of the data variables and methods of the parent class.
The child class can either utilize the methods as is or they can override the methods to modify their behavior without affecting the parent class or any objects that have instantiated that class.

Using inheritance in code development creates a hierarchy of objects, which often improves the readability of your code
It also saves time end effort by avoiding duplicate code production, i.e., inheriting from classes that have similar behavior and modifying them instead of writting a new class from scratch.

<!--- 
Take our previous Molecule class as an example.
We want to construct a new class for Molecules that prints with a different structure, however, assume we have code that utilized our original Molecule class, which means we cannot modify it.
We will inherit our original class and modify the print method to meet our new needs.
--->

### Interfaces and Abstract Classes
An interface is a way to define how a class will be designed without implementing any of the methods.
An interface can be inherited by another interface, in which case it is usually extended into a larger interface by the child interface, or it can be inherited by a class, which implements at least all of the methods from the interface.
An interface creates a structure for code to look for without defining how the methods are implemented.

Abstract Classes are interfaces that have implementations for one or more of their methods.
The implementations within an abstract class can be used by their children or overridden to achieve new behavior.

In python, we construct interfaces and abstract classes in a similar way.
We first need to import from the python library the abstract base class and abstract method decorator.
~~~
from abc import ABC, abstractmethod
~~~
{: .language-python}
`ABC`{: .language-python} is the abstract base class.
Inheriting from the abstract base class enforces that any child classes of the interface must have some implementation of any abstract methods.
~~~
class interface_sample(ABC):
~~~
{: .language-python}
Here we define the name of our interface.
Note `ABC`{: .language-python} in the class definition to denote the inheritance.
Then we create a definition for each method we would like any of the interfaces children to implement.
~~~
	@abstractmethod
	def first_method(self):
		pass
		
	@abstractmethod
	def second_method(self):
		pass
~~~
{: .language-python}
Since we have used the `@abstractmethod`{: .language-python} decorator for each of these methods, any child class must implement a method with the same name as the interface.

To convert this class from an Interface to an Abstract Class, we simply include one or more methods that are already implemented.
~~~
class abstract_class_sample(ABC):
	@abstractmethod
	def first_method(self):
		pass
		
	def second_method(self):
		print('Doing something in this method so it is implemented.')
~~~
{: .language-python}

<!---
-An Interface is a class definition that does not implement any of its methods. It is explicitly used as a contract to specify what methods any of its child classes will have to implement.

### Abstract Classes
-Abstract Classes are similar to interfaces. They define a base structure that subclasses inherit. The primary difference is that an interface has no method implementations, it simply has the class and method definitions. An abstract class can have implementations for one or more of its methods that can either be inherited or overridden by any of their subclasses.
--->

## Polymorphism
Polymorphism the concept of using different classes in place of one another.
Specifically, an object is polymorphic if it can be used in place of one or more classes or interfaces.
The intended use of a polymorphic object is to allow objects that are children of a parent class to be used as their parent class or for multiple objects that inherit from an interface to be used as the interface.

## Static Methods

For further reading on Object Oriented Programming in Python, consider this [tutorial](https://www.python-course.eu/python3_object_oriented_programming.php) or the Python3 documentation found [here](https://docs.python.org/3/tutorial/classes.html).
