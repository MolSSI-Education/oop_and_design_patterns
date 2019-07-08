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
molecule_name = "water molecule"
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

For this type of problem, we will want to create something called a *class*. Classes provide a way to bundle data and other functionality together. 

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
~~~
{: .language-python}

This is a simple definition of a class named Molecule. Let's look at what each line does.

The first line of this code
~~~
class Molecule:
~~~
{: .language-python}
is defining the name of the class as Molecule.
We then have a method called a constructor and it is called whenever you are instantiating an object of the class.
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
set the value of the local object variables to the value of the parameters. Here, the `self` syntax refers to the instance of the class. Any time you want to set or create a variable associated with a class in its definition, you use this syntax. 

We can now use this class definition in our code. For example, to create our water molecule, we use the class. This is called creating an *instance* of the class.

~~~
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"])
~~~
{: .language-python}

`mol1` in our code is now an object. We can access the variables associated with this instance of the molecule class using a dot notation. Variables associated with a class are also called `attributes`.

~~~
print(mol1.name)
print(mol1.charge)
print(mol1.symbols)
~~~
{: .language-python}

You should see the output

~~~
water dimer
0.0
['O', 'H', 'H']
~~~
{: .output}

> ## Check your understanding
> Create another instance of the class called `mol2`. This molecule should be an He molecule with 0 charge. After you have created this, print the molecule name and charge.
>> ## Solution
>> ~~~
>> mol2 = Molecule(name="He", charge=0.0, symbols=["He"])
>>
>> print(mol2.name)
>> print(mol2.charge)
>> ~~~
>> {: .language-python}
>>
>> ~~~
>> He
>> 0.0
>> ~~~
>> {: .language-output}
> {: .solution}
{: .challenge}

You may notice that if you print your molecules, you get something that is confusing and not so pretty.

~~~
print(mol1)
~~~
{: .language-python}

~~~
<__main__.Molecule object at 0x103f046d8>
~~~
{: .output}

We can create a nicer representation for printing by writing a `__str__` method for the class. In Python, there are special methods associated with classes which you can use for customization. These are also called "magic" methods. They exist inside a class, and begin and end with two underscores (`__`). The `__init__` we have already used is a magic method used to set initial properties of a class instance. The `__str__` method is called by built-in Python functions `print()` and `format()`. The return value of this function *must* be a string. 

The `__str__` method is simply a method to compute the string representation of our Molecule object to be used in printing, similar to how we defined it without any class, but it now will work for each instance of a Molecule without any modification. Let's add this to our class definition, making the whole definition look like the following

~~~
class Molecule:
    def __init__(self, name, charge, symbols):
        self.name = name
        self.charge = charge
        self.symbols = symbols
		
    def __str__(self):
        return 'name: ' + str(self.name) + '\ncharge:' + str(self.charge) + '\nsymbols:' + str(self.symbols)
		
mol1 = Molecule('water molecule', 0.0, ["O", "H", "H"])
mol2 = Molecule('He', 0.0, ["He"])
~~~
{: .language-python}

Now print these objects

~~~
print(mol1)
print(mol2)
~~~
{: .language-python}

You should see output that looks like this

~~~
name: water molecule
charge:0.0
symbols:['O', 'H', 'H']
name: He
charge:0.0
symbols:['He']
~~~
{: .language-output}

The construction and use of an object constructed through a class is simpler and more intuitive that trying to construct one from an arbitrary set of variables.
Now anytime we wish to create a new molecule object and print out its values, we only require two new lines of code.

> ## Check your understanding
> Add an additional attribute to the Molecule class (in `__init__`) which stores the number of atoms in the molecule. Print the number of atoms for the water molecule.
>> ## Solution
>> ~~~
>> class Molecule:
>>  def __init__(self, name, charge, symbols):
>>      self.name = name
>>      self.charge = charge
>>      self.symbols = symbols
>>      self.num_atoms = len(symbols)
>>
>>  def __str__(self):
>>      return 'name: ' + str(self.name) + '\ncharge:' + str(self.charge) + '\nsymbols:' + str(self.symbols)
>>
>> mol1 = Molecule('water molecule', 0.0, ["O", "H", "H"])
>>
>> print(F'{mol1.name} has {mol1.num_atoms} atoms.')
>> ~~~
>> {: .language-python}
>>
>> ~~~
>> water molecule has 3 atoms.
>> ~~~
>> {: .language-output}
> {: .solution}
{: .challenge}

<!--
## Composition

-including objects of class types within a class is called composition.
-We include instances of classes within
-->
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

Using inheritance in code development creates a hierarchy of objects, which often improves the readability of your code.
It also saves time end effort by avoiding duplicate code production, i.e., inheriting from classes that have similar behavior and modifying them instead of writting a new class from scratch.

Let us consider an example of a records system for a university.
A university has a large number of people, whether they are students or faculty.
We will start by creating some classes for each of these types of people.
First is a student class, at its simplest, a student has a name, a surname, and maybe a set of courses that they are registered for. Lets create a student class that takes a name and surname as parameters.
~~~
class Student:
    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.courses = []
~~~
{: .language-python}
Methods that only act upon the students data should be contained by the student class, this helps structure the methods in a more readible and accessible way.
Lets add three methods, a method to enroll a student in a course, a method to let the student drop a course, and the built in `__str__` method for the student class.
~~~
class Student:
    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.courses = []

    def enroll(self, new_course):
        self.courses.append(new_course)
~~~
{: .language-python}
It is often useful to generate a string representation of our class, so we want to override the built in method `__str__`.
~~~
class Student:
    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.courses = []

    def enroll(self, new_course):
        self.courses.append(new_course)
        
    def __str__(self):
        return self.surname + ", " + self.name + "\nCourses:\n" + str(courses)
~~~
{: .language-python}

> ## Check your understanding
> Add an additional method to the Student class to remove a course from the students enrolled courses.
>> ## Solution
>> ~~~
>> class Student:
>>   def __init__(self, name, surname):
>>      self.name = name
>>      self.surname = surname
>>      self.courses = []
>>
>>  def enroll(self, new_course):
>>      self.courses.append(new_course)
>>
>>  def drop_course(self, course):
>>      self.courses.remove(course)
>>
>>  def __str__(self):
>>      return self.surname + ", " + self.name + "\nCourses:\n" + str(courses)
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

Similar to the student, lets build a Faculty class to represent the instructors of the university.
Like the students, they have a name and a surname, but unlike the student they have a position denoting if they are a professor or lecturer and a salary.
~~~
class Faculty:
    def __init__(self, name, surname, position, salary):
        self.name = name
        self.surname = surname
        self.position = position
        self.salary = salary
        self.courses = []
        
    def __str__(self):
        return self.surname + ", " + self.name + "\nCourses:\n" + str(courses)
~~~
{: .language-python}
Like a student, a faculty has a set of courses, so we need to have methods to assign and unassign courses from their teaching load.
~~~
class Faculty:
    def __init__(self, name, surname, position, salary):
        self.name = name
        self.surname = surname
        self.position = position
        self.salary = salary
        self.courses = []
        
    def assign_course(self, new_course):
        self.courses.append(new_course)
    
    def unassign_course(self, course):
        self.courses.remove(course)
        
    def __str__(self):
        return self.surname + ", " + self.name + "\nCourses:\n" + str(courses)
~~~
{: .language-python}
Having built both a Student class and a Faculty class, notice the similarities between the two.
For variable, both classes have a name, a surname, and a set of courses.
For methods, both classes have a similar `__init__` method and a similar `__str__` method.
What if we want to add a new method to both classes? Consider a university ID number; most, if not all, universities generate id numbers for their students, faculty, and staff to avoid ambiguity that can arise from similar names.

If we want to add a new method to generate the id number of a given student or faculty, we have to add the method to both classes, which is duplicating the code in multiple places. This leads to more work for no tangible gain, not to mention, leads to multiple opportunities for mistakes to be made.
We can use inheritance to combat these problems.
We want to make a person class that contains the similarities of each class to act as their parent.
~~~
class Person:
    def __init__(self, name, surname):
        self.name = name
        self.surname = surname
        self.id = self.generate_id()
    
    def generate_id(self):
        id_hash = 0
        for s in self.name:
            id_hash += ord(s)
        for s in self.surname:
            id_hash = id_hash * ord(s)
        return id_hash % 1000000000
    
    def __str__(self):
        return self.surname + ", " + self.name + "\tID: " + str(self.id)
~~~
{: .language-python}

Now we can make the student class a child of the person class.
~~~
class Student(Person):
    def __init__(self, name, surname):
        self.courses = []
        super().__init__(name, surname)
    
    def __str__(self):
        return super().__str__() + "\nCourses:\n" + str(courses)
        
    def enroll(self, new_course):
        self.courses.append(new_course)
        
    def drop_course(self, course):
        self.courses.remove(course)
~~~
{: .language-python}
In both the `__init__` and `__str__` methods, we are using `super()`, which references the parent class of student, in this case Person, and calls the `__init__` and `__str__` methods to help initialize the class.
`super().__init__(name,surname)` tells python to call the `__init__` method of the parent class to initialize the name and surname variables. Now any changes to Person's `__init__` will also update the Student's `__init__`.

> ## Check your understanding
> Update the Faculty class to use the Person class as a parent.
>> ## Solution
>> ~~~
>> class Faculty(Person):
>>  def __init__(self, name, surname, position, salary):
>>      self.position = position
>>      self.salary = salary
>>      self.courses = []
>>      super().__init__(name, surname)
>>  
>>  def __str__(self):
>>      return super().__str__() + "\nCourses:\n" + self.courses
>>  
>>  def assign_course(self, new_course):
>>      self.courses.append(new_course)
>>  
>>  def unassign_course(self, course):
>>      self.courses.remove(course)
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

We initially created the person class to simplify the method to generate ids, but that method is not present in either class.
This is because both classes inherit the generate_id method from the Person class. Since they are not modifying the method, it does not need to appear within either child.
However, we can test the method to make sure it is working.
~~~
student1 = Student("John", "Smith")
print(student1)
~~~
{: .language-python}
gives the output:
~~~
Smith, John	ID: 546320160
Courses:
[]
~~~
{: .language-output}
The generate_id method is called in the `__init__` method of Person, so each Student and Faculty will autimatically generate their id upon initialization.

We can create further classes that inherit from person to cover different people at the university, such as Staff.

## Composition
TODO: Create a courses class and utilize it in Student and Faculty to hold more sophisticated course information/functionality.

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
