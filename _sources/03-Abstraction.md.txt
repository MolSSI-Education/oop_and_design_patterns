# Abstraction

````{admonition} Overview
:class: overview

Questions:
- What is Abstraction?
- Why should I care about Data Abstraction?

Objectives:
- Understand the concepts behind Data Abstraction.
````

## Abstraction

Abstraction is the concept of hiding implementation details from the user, allowing them to know how to use the code/class without knowing how it actually works or any implementation details.
For example,  when you use a Coffee machine, you interact with its interface, but you don't actually know how it is preparing the coffee inside.
Another example is that when a Web browser connects to the Internet, it interacts with the Operating system to get the connection, but it doesn't know if you are connecting using a dial-up, dsl, cable, etc.

There are many benefits of using abstraction:
1. Can have multiple implementations
2. Can build complex software by splitting functionality internally into steps, and only exposing one method to the user
3. Change implementation later without affecting the user by moving frequently changing code into separate methods.
4. Easier code collaboration since developers don't need to know the details of every class, only how to use it.
5. One of the main concepts that makes the code flexible and maintainable.

In Python, abstraction can be achieved by the use of private and public attributes and methods.

Generally, variables can be public, in which case they are directly modifiable, or private, meaning interaction with their values is only possible through internal class methods.
In python, there are no explicitly public or private variables, all variables are accessible within an object.
However, the predominantly accepted practice is to treate names prefixed with an underscore as non-public. (See [Python Private Variables](https://docs.python.org/3/tutorial/classes.html#private-variables))

### Private Data
Making data private is a way to protect a user from creating errors by hiding sensitive variables.
We can utilize our previously defined `Molecule` class to provide an example.
`Molecule` is currently defined as:
````{tab-set-code} 

```{code-block} python
class Molecule:
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
		
    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}'
```
````


We would like to provide users with a way to determine the number of atoms present in the molecule.
There are a few ways to accomplish this behavior.
First, we could add a `num_atoms` parameter to our `__init__` method and have the user pass it in.
Though it will work, it is slightly redundant as the number of atoms should be equivalent to the number of symbols the user has already passed in.

We could create a method that will calculate the number of atoms from the list of symbols and return it to the user.
In this instance, that will work well as it is a simple calculation and we can just re-generate the number as needed.

However, lets consider this as an example where re-generating every time it is needed is not a good idea, due to the complexity of the calculation.
In this case, we want to store the number of atoms so we can simply reference it as needed.

As a first attempy, why don't we just add a new variable to the `__init__` method that infers the value based on the list of symbols.
````{tab-set-code} 

```{code-block} python
class Molecule:
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates
        self.num_atoms = len(symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}\nnumber of atoms: {self.num_atoms}'

```
````

Recreating our two molecules using the new class and printing them, we get the following code and output:
````{tab-set-code} 

```{code-block} python
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
mol2 = Molecule(name="He", charge=0.0, symbols=["He"], coordinates=[[0,0,0]])
print(mol1)
print(mol2)
```
````

````{tab-set-code} 

```{code-block} output
name: water molecule
charge: 0.0
symbols: ['O', 'H', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
number of atoms: 3
name: He
charge: 0.0
symbols: ['He']
coordinates: [[0, 0, 0]]
number of atoms: 1
```
````

So far so good, but what if the list of symbols needs to change?
For the sake of the example, we will update the list of symbols to remove a hydrogen.
````{tab-set-code} 

```{code-block} python
mol1.symbols = ["O", "H"]
print(mol1)
```
````

````{tab-set-code} 

```{code-block} output
name: water molecule
charge: 0.0
symbols: ['O', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
number of atoms: 3
```
````


We can see that, though the list of symbols has properly updated, the number of atoms did not change.

## Property and Setter
We can utilize abstraction to cover this issue.
Python has two decorators, modifiers that can be applied to methods, that will abstract the behavior of the `num_atoms` variable: `property` and `setter`.
These allow attributes to be used in a pythonic way, while allowing more control over their values.

We want to update the `num_atoms` variable whenever the value of `symbols` is updated.
We can do this by creating a property and setter method for the 'symbols' variable.
````{tab-set-code} 

```{code-block} python
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
        self.num_atoms = len(symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}\nnumber of atoms: {self.num_atoms}'
```
````

To utilize the `property` decorator, we preempt a method with `@` followed by the decorator name.
The specific method used should have the same name as the variable we want the user to reference, in this case `symbols`.
Notice that the `symbols()` method is referencing a variable called `_symbols`, remember that prefacing a variable with an underscore infroms the user they should not edit this variable.

The `setter` decorator is applied in a similar way.
We use an `@` followed by the name of the referenced variable, in this case `symbols` and end with `.setter`.
This is preempting a method of the same name, `symbols`.

When a user tries to update the value of `symbols`, python will utilize the setter method to update its value, giving us more control over how that update occurs.
In this case, we are setting the value of `_symbols` to the value provided by the user, then updating the value of `num_atoms` based on the new value.
Note that we can remove the line that sets `num_atoms` from the constructor, as that will automatically occur when `symbols` is set.

If we recreate `mol1` and attempt to edit its `symbols` now we should see a different result:
````{tab-set-code} 

```{code-block} python
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
mol1.symbols = ["O", "H"]
print(mol1)
```
````

````{tab-set-code} 

```{code-block} output
name: water molecule
charge: 0.0
symbols: ['O', 'H']
coordinates: [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
number of atoms: 2
```
````


Now, whenever someone tries to access 'symbols', it will properly return the value stored in our private variable.
When someone tries to update the value of symbols directly, it will update the private variable, and also update the 'num_atoms' variable.

### Private Methods
Like variables, it can be useful to have private methods, methods the user should not directly call.
Often these can be helper methods which are performing portions of a calculation that are not useful when run outside of the wider calculation.
To include private methods in a class, we similarly preempt the method name with the `_` character. This does not explicitly make the function private, but by convention informs other developers and users that this method should not be called directly.

For the sake of showing a simple example, let us assume that updating the value of `num_atoms` is a much more complicated procedure than it is under the current definition of `Molecule`.
Since we want to keep the content of the `symbols` setter method fairly straightforward and easy to understand, we create a helper method to handle the updating of `num_atoms`.
````{tab-set-code} 

```{code-block} python
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
```
````

Now when a user updates the value of `symbols`, python will call the *private* method `_update_num_atoms()`, which will correctly update the number of atoms stored in the object.
In this instance it is not necessary to hide the method `_update_num_atoms()` as no harm will occur if it is directly called, but it works as a useful example of how a method can be hidden.

## Name Mangling
For the previous sections, we have utilized python conventions to make our variables private.
By prefacing a variable or method with and underscore, `_`, we are telling any users not to directly interact with the variable or method.
However, if we look at the dictionary of the object, we can see that they are still accessible and users may try and use them.
````{tab-set-code} 

```{code-block} python
print(dir(mol1))
```
````

````{tab-set-code} 

```{code-block} output
['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_symbols', '_update_num_atoms', 'charge', 'coordinates', 'name', 'num_atoms', 'symbols']
```
````


Python allows us to further hide them through a process called name mangling.
To apply name mangling to a variable or method, simply preface a value with two underscores `__`.
Let us try and apply name mangling to both the `_symbols` variable and the `_update_num_atoms` method.
````{tab-set-code} 

```{code-block} python
class Molecule:
    def __init__(self, name, charge, symbols, coordinates):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.coordinates = coordinates

    @property
    def symbols(self):
        return self.__symbols
        
    @symbols.setter
    def symbols(self, symbols):
        self.__symbols = symbols
        self.__update_num_atoms()

    def __update_num_atoms(self):
        self.num_atoms = len(self.symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}\ncoordinates: {self.coordinates}\nnumber of atoms: {self.num_atoms}'
```
````

We can then recreate `mol1` and check its dictionary.
````{tab-set-code} 

```{code-block} python
mol1 = Molecule(name='water molecule', charge=0.0, symbols=["O", "H", "H"], coordinates=[[0,0,0],[0,1,0],[0,0,1]])
print(dir(mol1))
```
````

````{tab-set-code} 

```{code-block} output
['_Molecule__symbols', '_Molecule__update_num_atoms', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'charge', 'coordinates', 'name', 'num_atoms', 'symbols']
```
````

We can now see that both the private variable and method have been renamed internally and prefaced with their class name: `_Molecule__symbols` and `_Molecule__update_num_atoms`.
Though they are still visible, this further pushes them behind a layer of obfuscation, strongly suggesting to any user that they should not directly interact with these values.
This is as close to making data private as python will currently allow.
````{admonition} Key Points
:class: key

- Public vs Private variables and methods
- Python syntax for private variables and methods
````
