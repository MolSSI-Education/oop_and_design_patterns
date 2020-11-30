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

We will use the three classes defined in the `Inheritance` lessons, `Person`, `Faculty`, and `Student` to create an example of this behaviour.

First we will restate the three classes here. First the parent class Person.

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
            id_hash *= ord(s)
        return id_hash % 1000000000
    
    def __str__(self):
        return f'{self.surname}, {self.name}\tID: {self.id}'
~~~
{: .language-python}

Then the two child classes, Faculty and Student.
~~~
class Faculty(Person):
    def __init__(self, name, surname, position, salary):
        self.position = position
        self.salary = salary
        self.courses = []
        super().__init__(name, surname)

    def __str__(self):
        return super().__str__() + f'\nCourses:\n{self.courses}'

    def assign_course(self, new_course):
        self.courses.append(new_course)

    def unassign_course(self, course):
        self.courses.remove(course)

class Student(Person):
    def __init__(self, name, surname):
        self.courses = []
        super().__init__(name, surname)
    
    def __str__(self):
        return super().__str__() + f'\nCourses:\n{self.courses}'
        
    def enroll(self, new_course):
        self.courses.append(new_course)
        
    def drop_course(self, course):
        self.courses.remove(course)
~~~
{: .language-python}

Since `Faculty` and `Student` are both children of the `Person` class, take note that they both share a `name`, `surname`, and a `id` with their parent.
To ensure that both children are polymorphic, we want to ensure that any method that operates on a `Person` will correctly operate on them as well.

Here we will build a simple method that utilizes variables of a `Person` to provide a formatted output.

~~~
def formatted_person_print(person):
    print(f'{person.id}\n{person.name} {person.surname}')
~~~
{: .language-python}

We will create a person use them in our method.

~~~
person_example = Person('John', 'Smith')
print(person_example)
~~~
{: .language-python}
~~~
Smith, John     ID: 546320160
~~~
{: .language-output}

Using out new method on this person we get our new format.
~~~
formatted_person_print(person_example)
~~~
{: .language-python}
~~~
546320160
John Smith
~~~
{: .language-output}

For `Faculty` and `Student` to be polymorphic, they should both work with this method.

~~~
faculty_example = Faculty('Jane', 'Doe', 'Department Chair', 160000)
student_example = Student('James', 'Smith')
formatted_person_print(faculty_example)
formatted_person_print(student_example)
~~~
{: .language-python}
~~~
291216936
Jane Doe
167856640
James Smith
~~~
{: .language-output}

We have properly made the two children polymorphic. Any method that utilizes `Person` objects should be able to use objects of either `Faculty` or `Student`. This allows us to extend the behaviour of a `Person` without breaking any code that relies on it.