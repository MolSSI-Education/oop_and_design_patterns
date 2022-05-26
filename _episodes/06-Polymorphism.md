---
title: "Polymorphism"
teaching: 0
exercises: 0
questions:
- "What is Polymorphism?"
objectives:
- "Understand the concepts behind Polymorphism."
keypoints:
- "Encapsulation"
- "Inheritance"
---

{% include links.md %}
## Polymorphism
Polymorphism the concept of using different classes in place of one another.
Specifically, an object is polymorphic if it can be used in place of one or more classes or interfaces.
The intended use of a polymorphic object is to allow objects that are children of a parent class to be used as their parent class or for multiple objects that inherit from an interface to be used as the interface.

This is very similar to python's concept of `Duck Typing`, i.e. if it looks like a duck and quacks like a duck, then it must be a duck.
Duck typing is more of a passive concept, we assume objects of certain types can be used in certain ways and simply try to use them. For example, if we believe a variable is a number, we assume we can perform mathematical operations with it.

Polymorphism is the practice of making sure duck typing will work.

We will use the two classes, `Molecule` and `AtomMolecule` to create an example of polymorphism.

First we will restate the two classes here.
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

Since `AtomMolecule` is a child of the `Molecule` class, take note that they both share a many variable names..
To ensure that `AtomMolecule` is polymorphic, we want to ensure that any method that operates on a `Molecule` will correctly operate on an instance of `AtomMolecule` as well.

Here we will build a simple method that utilizes variables of a `Molecule` to provide a formatted output.

~~~
def formatted_print(molecule):
    return f'{molecule.name} is made of {molecule.symbols} and has an atomic charge of {molecule.charge}'
~~~
{: .language-python}

We will create a `Molecule` to provide to the method.

~~~
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
formatted_print(mol1)
~~~
{: .language-python}
~~~
"water molecule is made of ['O', 'H', 'H'] and has an atomic charge of 0.0"
~~~
{: .language-output}

For `AtomMolecule` to be polymorphic, it should also work with the method.

~~~
oxygen = Atom("oxygen", "O", 8, 15.999, [0,0,0])
hydrogen1 = Atom("hydrogen", "H", 1, 1.00784, [0,1,0])
hydrogen2 = Atom("hydrogen", "H", 1, 1.00784, [0,0,1])

mol2 = AtomMolecule(name='water molecule', charge=0.0, atoms=[oxygen, hydrogen1, hydrogen2])
formatted_print(mol2)
~~~
{: .language-python}
~~~
"water molecule is made of ['O', 'H', 'H'] and has an atomic charge of 0.0"
~~~
{: .language-output}

We have properly made `AtomMolecule` polymorphic. Any method that utilizes `Molecule` objects should be able to use `AtomMolecule` objects. This allows us to extend the behaviour of a `Molecule` without breaking any code that relies on it.