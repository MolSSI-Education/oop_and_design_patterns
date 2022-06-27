---
title: "Inheritance"
teaching: 0
exercises: 0
questions:
- "What is Inheritance?"
- "Why should I use Inheritance in my code?"
objectives:
- "Understand the concepts behind Inheritance."
keypoints:
- "Parent vs Child classes"
- ""
---

## Inheritance
Inheritance is the principle of extending a class to add capabilities without modifying the original class.
We call the class that is being inherited the parent, and the class that is inheriting the child.
The child class obtains the properties and behaviors of its parent unless it overrides them.

In coding terms, this means a class that inherits from a parent class will contain all of the data variables and methods of the parent class by default.
The child class can either utilize the methods as is or they can override the methods to modify their behavior without affecting the parent class or any objects that have instantiated that class.

Using inheritance in code development creates a hierarchy of objects, which often improves the readability of your code.
It also saves time end effort by avoiding duplicate code production, i.e., inheriting from classes that have similar behavior and modifying them instead of writting a new class from scratch.

We can utilize our `Molecule` class to create an example.
Of note, we are using the version prior to applying name mangling.
Since name mangling is applied based on the class it is defined in and not where it is called from, it will cause issues with inheritance.
~~~
class Molecule:
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
~~~
{: .language-python}

Currently the `Molecule` class takes in lists for both `symbols` and `coordinates`, relying on the indexes of each list to tie the appropriate atom symbol to it's coordinates.
Depending on your software's needs, it may be more appropriate to store each atom as an object, allowing the `Molecule` class to store a list of atoms.

A simple `Atom` class may look like the following:
~~~
class Atom:
    def __init__(self, name, symbol, number, mass, coordinates):
        self.name = name
        self.symbol = symbol
        self.number = number
        self.mass = mass
        self.coordinates = coordinates
~~~
{: .language-python}

Given this `Atom` class, we have a few possibilities for incorporating it into the Molecule class.
We could try and make the `symbols`, `coordinates`, and `atoms` variables optional and rely on keyword arguments to properly assign them in the constructor.
Alternatively we could utilize inheritance to extend the behavior of the `Molecule` class while preserving its behavior.

To do so, we can create a new class called `AtomMolecule`.
First, let us look at just the class definition and the constructor.
~~~
class AtomMolecule(Molecule):
    def __init__(self, name, charge, atoms):
        self.atoms = atoms
        super().__init__(name, charge, self.symbols, self.coordinates)
~~~
{: .language-python}

Notice that after the class name, we have now included `(Molecule)`.
Python uses the above syntax to specify inheritance.
The class definition specifies that `AtomMolecule` should inherit the data and methods from its parent class `Molecule`.

As far as the constructor is concerned, we are not passing in a set of symbols or coordinates, as those are both derived from the set of atoms.
We set the given list of atoms to the instance variable `self.atoms`.
Finally, we are utilizing the `__init__()` method from the parent class.
The syntax `super().` is telling python where to look to find `__init__()`, so that it searches in the parent class.

If you try and run just this code and create an `AtomMolecule`, you will run into some errors.
Even though we are inheriting from `Molecule`, `AtomMolecule` has no understanding of `symbols`, so we cannot use the constructor.
A solution to this is to add a setter and property for atoms so we can control how `atoms` is updated.
When `atoms` is updated, we need to update both `symbols` and `coordinates`.
Here is one option for how to generate the values for `symbols`:
~~~
def _update_symbols(self):
        list_symbols = []
        for atom in self.atoms:
            list_symbols.append(atom.symbol)
        self._symbols = list_symbols
~~~
{: .language-python}
which iterates through each atom and appends the symbol to the list.
Adding a similar function for `coordinates` creates the complete class.
~~~
class AtomMolecule(Molecule):
    def __init__(self, name, charge, atoms):
        self.atoms = atoms
        super().__init__(name, charge, self.symbols, self.coordinates)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._atoms = atoms
        self._update_symbols()
        self._update_coordinates()

    def _update_symbols(self):
        list_symbols = []
        for atom in self.atoms:
            list_symbols.append(atom.symbol)
        self._symbols = list_symbols

    def _update_coordinates(self):
        list_coordinates = []
        for atom in self.atoms:
            list_coordinates.append(atom.coordinates)
        self.coordinates = list_coordinates
~~~
{: .language-python}

To test it we can create a set of atoms to pass into an `AtomMolecule`:
~~~
oxygen = Atom("oxygen", "O", 8, 15.999, [0,0,0])
hydrogen1 = Atom("hydrogen", "H", 1, 1.00784, [0,1,0])
hydrogen2 = Atom("hydrogen", "H", 1, 1.00784, [0,0,1])

mol1 = AtomMolecule(name='water molecule', charge=0.0, atoms=[oxygen, hydrogen1, hydrogen2])
print(mol1)
~~~
{: .language-python}
~~~
name: water molecule
charge: 0.0
symbols: ['O', 'H', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
number of atoms: 3
~~~
{: .language-output}

We can see in the output that `AtomMolecule` is correctly printing out the contents of the molecule, but we have not defined a `__str__()` method.
Since `AtomMolecule` is inheriting from `Molecule`, it is inheriting all of the methods defined in `Molecule`.
As a dictionary based language, when you call a method, python will first look in the object you are calling the method from.
If it is unable to find the method, it will attempt to look in the parent class of the object.
It will continue looking up through the inheritance tree until it either finds the method or runs out of parents to look in.

Similarly, we did not add a setter for `symbols` in `AtomMolecule`, yet the number of atoms was correctly set, since we inherited the property and setter methods from `Molecule`.

By creating `AtomMolecule`, we have extended the behavior of `Molecule`.
Instead of only working with a set of symbols and coordinates, we are now capable of working with a set of `Atom` objects.
Note that the external behavior/functionality of `AtomMolecule` is still identical to `Molecule`, but the internal operations are different to adjust to different input.

### Composition and Aggregation
When creating `AtomMolecule`, we used a type of object inside another object.
We grouped a set of related data together through encapsulation, then utilized the created object in another class.
There are two different forms that this can take, Composition and Aggregation.
The main difference is the ownership of the object.

With the current example of `AtomMolecule`, we are creating `Atoms` and passing them into the object.
If the `AtomMolecule` object was deleted, the defined `Atoms` would persist.
This is an example of aggregation.
The `AtomMolecule` object is using the `Atoms`, but it does not have ownership of them.

If, instead, we had the `Molecule` class create atoms from the set of symbols and coordinates, the `Molecule` class would own the `Atoms`, which would be an example of composition.