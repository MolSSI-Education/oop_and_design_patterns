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
Encapsulation is the concept of enclosing related data and methods acting on those data within a single unit called a class.
A class will consist of a set of data (variables) and a set of methods that interact with the data.
There are a number of benefits to creating classes.
1. It aids in the understanding of the code being developed.
It is much easier to understand how an object will behave if all the data and methods that interact with that data are enclosed within its class.
It means all of the information about that object is grouped into a single location.
2. As an extension of the first benefit, anyone developing code that utilizes your objects will have a clear understanding of how they are allowed to interact with it.
The methods they can use will be located within the class of the object.
3. It promotes security of data by restricting the ways to access data to the specific methods within the class.
4. Methods have full access to their data so that you don’t have to keep passing data and parameters between methods. Also, this way, you avoid the use of global variables.


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

print(f'name: {molecule_name}\ncharge: {molecule_charge}\nsymbols: {molecule_symbols}')

molecule2_name = "He"
molecule2_charge = 0.0
molecule2_symbols = ["He"]

print(f'name: {molecule2_name}\ncharge: {molecule2_charge}\nsymbols: {molecule2_symbols}')

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
water molecule
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
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}'
		
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
charge: 0.0
symbols: ['O', 'H', 'H']
name: He
charge: 0.0
symbols: ['He']
~~~
{: .language-output}

The construction and use of an object constructed through a class is simpler and more intuitive that trying to construct one from an arbitrary set of variables.
Now anytime we wish to create a new molecule object and print out its values, we only require two new lines of code.

> ## Check your understanding
> Add an additional attribute to the Molecule class (in `__init__`) which stores the number of atoms in the molecule. Print the number of atoms for the water molecule.
>> ## Solution
>> ~~~
>> class Molecule:
>>     def __init__(self, name, charge, symbols):
>>         self.name = name
>>         self.charge = charge
>>         self.symbols = symbols
>>         self.num_atoms = len(symbols)
>>
>>     def __str__(self):
>>         return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}'
>> 
>> mol1 = Molecule('water molecule', 0.0, ["O", "H", "H"])
>>
>> print(f'{mol1.name} has {mol1.num_atoms} atoms.')
>> ~~~
>> {: .language-python}
>>
>> ~~~
>> water molecule has 3 atoms.
>> ~~~
>> {: .language-output}
> {: .solution}
{: .challenge}