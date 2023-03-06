---
title: "Data Abstraction"
teaching: 0
exercises: 0
questions:
- "What is Data Abstraction?"
- "Why should I care about Data Abstraction?"
objectives:
- "Understand the concepts behind Data Abstraction."
keypoints:
- "Public vs Private variables and methods"
- "Python syntax for private variables and methods"
---

{% include links.md %}
## Data Abstraction

<!--- Private vs Public --->
<!--Data abstraction is the concept of hiding variables behind methods to simplify code, increase development efficiency and improve security. -->

Data abstraction is the concept of hiding implementation details from the user, allowing them to know how to use the code/class without knowing how it actually works or any implementation details. For example,  when you use a Coffee machine, you interact with its interface, but you don't actually know how it is preparing the coffee inside. Another example is that when a Web browser connects to the Internet, it interacts with the Operating system to get the connection, but it doesn't know if you are connecting using a dial-up, dsl, cable, etc.

Clearly there are many benefits of using data abstraction:
1. Can have multiple implementations
2. Can build complex software by splitting functionality internally into steps, and only exposing one method to the user
3. Change implementation later without affecting the user by moving frequently changing code into separate methods.
4. Easier code collaboration since developers don't need to know the details of every class, only how to use it.
5. One of the main concepts that makes the code flexible and maintainable.

In Python, data abstraction can be achieved by the use of private and public attributes and methods.


Generally, variables can be public, in which case they are directly modifiable, or private, meaning interaction with their values is only possible through internal class methods.
In python, there are no explicitly public or private variables, all variables are accessible within an object.
However, the predominantly accepted practice is to treate names prefixed with an underscore as non-public. (See [Python Private Variables](https://docs.python.org/3/tutorial/classes.html#private-variables))

### Private Data
The best way to achieve data abstraction of a variable in python is to use the '@property' decorator, and the 'setter' decorator. These allow attributes to be used in a pythonic way, while allowing more control over their values.
Consider our Molecule class:
~~~
class Molecule:
    def __init__(self, name, charge, symbols):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.num_atoms = len(symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}'
~~~
{: .language-python}

We already can see a prime candidate for this approach in our 'symbols' and 'num_atoms'. Since 'num_atoms' is based on the number of symbols, we dont want the user to modify it, and in addition, we want to update it whenever the number of symbols is updated. We can do this by creating a property and setter method for the 'symbols' variable.
~~~
class Molecule:
    def __init__(self, name, charge, symbols):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.num_atoms = len(symbols)

    @property
    def symbols(self):
        return self._symbols
        
    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols
        self.num_atoms = len(symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}'
~~~
{: .language-python}

Now, whenever someone tries to access 'symbols', it will properly return the value stored in our private variable. When someone tries to update the value of symbols directly, it will update the private variable, and also update the 'num_atoms' variable.

### Private Methods
To include private methods in a class, we similarly preceed the variable name with the `_` character. This does not explicitly make the function private, but by convention informs other developers and users that this method should not be called directly.

Consider we are writting an algorithm for a computational chemistry method and we wish to separate it into multiple stages for ease of use. We have stages of the algorithm that, if called by an external user at the wrong time, could corrupt the accuracy of out data objects and cause incorrect results.

We would like to hide the internal methods, and only allow the method call for the algorithm itself to be used by the user.

~~~
class Molecule:
    def __init__(self, name, charge, symbols):
        self.name = name
        self.charge = charge
        self.symbols = symbols
        self.num_atoms = len(symbols)

    @property
    def symbols(self):
        return self._symbols
        
    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols
        self.num_atoms = len(symbols)

    def __str__(self):
        return f'name: {self.name}\ncharge: {self.charge}\nsymbols: {self.symbols}'
        
    def my_algorithm(self, parameters_list):
        # Performs our algorithm.
        _first_stage(params)
        _second_stage(params)
        _third_stage(params)
        
    def _first_stage(self, parameters_list):
        # Performs the first stage of the algorithm.
    
    def _second_stage(self, parameters_list):
        # Performs the second stage of the algorithm.
    
    def _third_stage(self, parameters_list):
        # Performs the third stage of the algorithm.
~~~
{: .language-python}

We have separated out each of the stages out into its own function, prefacing their name with the `_` character to let users and developers know to either not call this method directly, or to use it very carefully. This is the best we can do to hide the methods within this class in python.
We have left out the `_` character in the method `my_algorithm` to let users know that it can be called by anyone.

If we need to change the implementation of the algorithm later down the road, we can modify any of the stages, or even add more stages, without changing how the algorithm is called. This is a great way to update important methods without breaking any code that relies upon them for their operation.