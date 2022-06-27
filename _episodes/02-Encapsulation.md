---
title: "Encapsulation"
teaching: 0
exercises: 0
questions:
- "What is encapsulation?"
- "Why should I use encapsulation in my code?"
objectives:
- "Understand the concept of a Class."
- "Understand the python syntax for Classes."
keypoints:
- "Encapsulation is the practice of enclosing data and functions in a class."
---

{% include links.md %}
## Encapsulation
Encapsulation is the concept of enclosing related data and methods acting on those data within a single object called a class.
A class will consist of a set of data (variables) and a set of methods that interact with the data.
There are a number of benefits to creating classes.
1. It aids in the understanding of the code being developed.
It is much easier to understand how an object will behave if all the data and methods that interact with that data are enclosed within its class.
It means all of the information about that object is grouped into a single location.
2. As an extension of the first benefit, anyone developing code that utilizes your objects will have a clear understanding of how they are allowed to interact with it.
The methods they can use will be located within the class of the object.
3. It enables data security.
A developer can restrict access to data within a class through specific methods.
4. Methods have full access to their data so that you don’t have to keep passing data and parameters between methods.
Also, this way, you avoid the use of global variables, which while useful, can lead to some issues such as name collision.


Classes are used in code to provide a general definition of an object.
We call an object an instance of a class, meaning it takes the structure of the class and fills in any necessary data.
Using classes as a general definition allows you to construct many instances of the class with little work.
Consider the case of creating a molecule.
First we decide on a few simple peices of data to make up a molecule:
* name - the name of our molecule
* charge - the molecular charge of the molecule.
* symbols - the list of atomic symbols that make up the molecule
* coordinates - the positions of each of the atoms.

Here we create a water molecule using just python variables:
~~~
name = "water molecule"
charge = 0.0
symbols = ["O", "H", "H"]
coordinates = [[0,0,0],[0,1,0],[0,0,1]]
~~~
{: .language-python}

Notice the first issue, we now have a set of four variables that, if seen separately, would not have nay connection to one another.
As the developer, we know that the symbols and coordinates are tied to the same molecule, but to any observer they are unrelated variables.

We can try and mitigate this through variable names:
~~~
molecule_name = "water molecule"
molecule_charge = 0.0
molecule_symbols = ["O", "H", "H"]
molecule_coords = [[0,0,0],[0,1,0],[0,0,1]]
~~~
{: .language-python}

Now at a glance they are all tied to a molecule.
We can add a little utility to our molecule by printing it out.
~~~
print(f'name: {molecule_name}\ncharge: {molecule_charge}\nsymbols: {molecule_symbols}\ncoordinates: {molecule_coords}')
~~~
{: .language-python}

Having created a single molecule, we can start to think of issues with this design.
It is likely that our code will need to deal with more than one molecule.

~~~
molecule2_name = "He"
molecule2_charge = 0.0
molecule2_symbols = ["He"]
molecule2_coords = [[0,0,0]]

print(f'name: {molecule2_name}\ncharge: {molecule2_charge}\nsymbols: {molecule2_symbols}\ncoordinates: {molecule2_coords}')
~~~
{: .language-python}

Under the current design, we need to repeat all of our previous code to create each additional molecule.
For a single instance, this may not seem terrible, but what if we need to make tens or hundreds of instances of the same type of object?
The amount of additional code that needs to be written grows very quickly.
In addition, the room for errors with naming grows for each new instance.
Each molecule has four variables that need to be correctly named or they will not make sense.

Encapsulation is designed to solve this problem.
We can create something called a *class* to hold our molecule information.
Classes provide a way to bundle data and other functionality together. 

## Class Syntax
Python allows for the creation of classes. We define a `Molecule` class.
~~~
class Molecule:
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
~~~
{: .language-python}

`Molecule` is a class that contains all the pieces of data we are currently associating with molecules.
Let's look at what each line does.

The first line of this code
~~~
class Molecule:
~~~
{: .language-python}
is defining the name of the class as Molecule.
We then have a method called a constructor and it is called whenever you are instantiating an object of the class.
We have the definition of the constructor
~~~
    def __init__(self, name, charge, symbols, coordinates):
~~~
{: .language-python}
that has four parameters, one for each piece of data that makes up a molecule: name, charge,  symbols, and coordinates.
The 5th parameter, `self` refers to the instance of the class.
Every method of a class must have a reference to the instance as the first variable.
This variable can be given any name, but by convention is usually called `self`.
The parameters of a constructor are required anytime you want to create an instance of the class. 
These can have default values if they are non-required variables.
The next four lines
~~~
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
~~~
{: .language-python}
set the value of the local object variables to the value of the parameters passed into the constructor.
Any time you want to set or create a variable associated with a class in its definition, you use this syntax. 
We can now use this class definition in our code.
For example, to create our water molecule, we use the class.
This is called creating an *instance* of the class.

~~~
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
~~~
{: .language-python}

`mol1` in our code is now an object.
Note that we do not pass a value to be assigned to `self` as python will automatically fill in the value for that variable when it calls `__init__()`.
We can access the variables associated with this instance of the molecule class using a dot notation.
Variables associated with a class are also called `attributes`.

~~~
print(mol1.charge)
print(mol1.symbols)
print(mol1.coordinates)
~~~
{: .language-python}

You should see the output

~~~
water molecule
0.0
['O', 'H', 'H']
[[0, 0, 0], [0, 1, 0], [0, 0, 1]]
~~~
{: .output}

> ## Check your understanding
> Create another instance of the class called `mol2`.
This molecule should be a recreation of the above He molecule.
After you have created this, print the molecule name and charge.
>> ## Solution
>> ~~~
>> mol2 = Molecule(name="He", charge=0.0, symbols=["He"], coordinates=[[0,0,0]])
>> ~~~
>> {: .language-python}
>> Following the syntax from creating the first molecule, we assign a variable name, and set it equal to a call to the Molecule class.
We then assign each value to the variables in the constructor.
We can then simply print both values.
>>
>> ~~~
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

## Printing Objects

You may notice that if you print your molecules, you get something that is confusing and not so pretty.

~~~
print(mol1)
~~~
{: .language-python}

~~~
<__main__.Molecule object at 0x7f6dec721f10>
~~~
{: .output}

The specific output may be different but the structure should be the same.

Python is printing the location of the object named `mol1`.
We can create a nicer representation for printing by writing a `__str__` method for the class.
In Python, there are special methods associated with classes which you can use for customization.
These are also called *magic* methods.
They exist inside a class, and begin and end with two underscores (`__`).
The `__init__` we have already used is a magic method used to set initial properties of a class instance.
The `__str__` method is called by built-in Python functions `print()` and `format()`. The return value of this function *must* be a string. 

The `__str__` method is simply a method to generate the string representation of our Molecule object to be used in printing, similar to how we defined it without any class, but it now will work for each instance of a Molecule without any modification.
Let's add this to our class definition, making the whole definition look like the following

~~~
class Molecule:
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
		
    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}'
		
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
mol2 = Molecule(name="He", charge=0.0, symbols=["He"], coordinates=[[0,0,0]])
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
charge: 0.0
symbols: ['O', 'H', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
name: He
charge: 0.0
symbols: ['He']
coordinates: [[0, 0, 0]]
~~~
{: .language-output}

Creating a new molecule object now requires only one line of code, with only one new variable name created and assigned.
This removes much of the redundancy of creating multiple variables and removes many possible points where a syntax error could occur.

In summary, utilizing encapsulation to wrap up the data and methods that act on them has provided a number of benefits.
We have reduced developer work by removing redundancy.
We have reduced the likelihood of errors appearing in the syntax written by reducing the number of possible locations for human error.
The next lesson will cover abstraction to show how encapsulation can improve the security of your code.