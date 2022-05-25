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
## Abstract Classes
An abstract class is a way to define how a class will be designed while allowing some or all of its methods to remain unimplemented.
An abstract class can be inherited by another abstract class, in which case it is usually extended into a larger abstract class by the child abstract class, or it can be inherited by a class, which implements at least all of the methods from the abstract class.
An abstract class creates a structure for code to look for without defining how the methods are implemented.
Any methods implemented within an abstract class can be used by their children or overridden to achieve new behavior.

In python, we construct abstract classes using a package from the standard library.
We first need to import from the python library the abstract base class and abstract method decorator.
~~~
from abc import ABC, abstractmethod
~~~
{: .language-python}
`ABC`{: .language-python} is the abstract base class.
Inheriting from the abstract base class enforces that any child classes of the abstract class must have some implementation of any abstract methods.
~~~
class AbstractClass(ABC):
~~~
{: .language-python}
Here we define the name of our abstract class.
Note `ABC`{: .language-python} in the class definition to denote the inheritance.
Then we create a definition for each method we would like any of the abstract classes children to implement and any methods we would like to pass on.
~~~
    @abstractmethod
    def first_method(self):
        pass

    @abstractmethod
    def second_method(self):
        pass
        
    def third_method(self):
        print('This method is implemented for use by any child classes.')
~~~
{: .language-python}
Since we have used the `@abstractmethod`{: .language-python} decorator for some of these methods, any child class must implement a method with the same name as the abstract class.

If we try and create an instance of this object, we will notice a problem.
~~~
abstract_class = AbstractClass()
~~~
{: .language-python}
~~~
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: Can't instantiate abstract class AbstractClass with abstract methods first_method, second_method
~~~~
{: .language-output}
Since there are unimplemented methods, we canot create the instance of the object.

A child class that inherits our abstract class will need to override `first_method` and `second_method` if it is not itself an abstract class.

~~~
class ChildClass(AbstractClass):
    def first_method(self):
        print('Implementing the first method inherited from my parent.')
    
    def second_method(self):
        print('Implementing the second method inherited from my parent.')
~~~
{: .language-python}

Note that we did not implement the `third_method` in the child class. Since `third_method` is not an abstract method, we can utilize the implementation provided by our parent.

Lets try and create and instance of the child class.
~~~
child_class = ChildClass()
~~~
{: .language-python}

We should see a blank output as the object has been created correctly.
We can now use all three methods on our child object.
~~~
child_class.first_method()
child_class.second_method()
child_class.third_method()
~~~
{: .language-python}
~~~
Implementing the first method inherited from my parent.
Implementing the second method inherited from my parent.
This method is implemented for use by any child classes.
~~~
{: .language-output}

If we wish to override the method provided by the abstract class, we can simply re-write the method.
~~~
class ChildClass(AbstractClass):
    def first_method(self):
        print('Implementing the first method inherited from my parent.')
    
    def second_method(self):
        print('Implementing the second method inherited from my parent.')
    
    def third_method(self):
        print('Overriding the third method provided by my parent.')
~~~
{: .language-python}

Re-instantiating the object and running the three methods, we will see the overridden output.
~~~
child_class = ChildClass()

child_class.first_method()
child_class.second_method()
child_class.third_method()
~~~
{: .language-python}
~~~
Implementing the first method inherited from my parent.
Implementing the second method inherited from my parent.
Overriding the third method provided by my parent.
~~~
{: .language-output}