---
title: "Abstract Classes"
teaching: 0
exercises: 0
questions:
- "What are Abstract Classes?"
- "Why should I use Abstract Classes?"
objectives:
- "Understand the concepts behind Abstract Classes."
keypoints:
- "Abstract Classes"
- "Abstract Methods"
- "Overwritting Abstract Methods"
---

{% include links.md %}
## Interfaces and Abstract Classes
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