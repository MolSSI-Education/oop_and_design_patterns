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
        return f'{self.surname}, {self.name}\nCourses:\n{self.courses}'
~~~
{: .language-python}

> ## Check your understanding
> Add an additional method to the Student class to remove a course from the students enrolled courses.
>> ## Solution
>> ~~~
>> class Student:
>>     def __init__(self, name, surname):
>>         self.name = name
>>         self.surname = surname
>>         self.courses = []
>>
>>     def enroll(self, new_course):
>>         self.courses.append(new_course)
>>
>>     def drop_course(self, course):
>>         self.courses.remove(course)
>>
>>     def __str__(self):
>>         return f'{self.surname}, {self.name}\nCourses:\n{self.courses}'
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
        return f'{self.surname}, {self.name}\nCourses:\n{self.courses}'
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
        return f'{self.surname}, {self.name}\nCourses:\n{self.courses}'
~~~
{: .language-python}
Having built both a Student class and a Faculty class, notice the similarities between the two.
For variables, both classes have a name, a surname, and a set of courses.
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
       return hash(self.name+self.surname) % 1000000000
    
    def __str__(self):
        return f'{self.surname}, {self.name}\tID: {self.id}'
~~~
{: .language-python}

Now we can make the student class a child of the person class.
~~~
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
In both the `__init__` and `__str__` methods, we are using `super()`, which references the parent class of student, in this case Person, and calls the `__init__` and `__str__` methods to help initialize the class.
`super().__init__(name,surname)` tells python to call the `__init__` method of the parent class to initialize the name and surname variables. Now any changes to Person's `__init__` will also update the Student's `__init__`.

> ## Check your understanding
> Update the Faculty class to use the Person class as a parent.
>> ## Solution
>> ~~~
>> class Faculty(Person):
>>     def __init__(self, name, surname, position, salary):
>>         self.position = position
>>         self.salary = salary
>>         self.courses = []
>>         super().__init__(name, surname)
>>  
>>     def __str__(self):
>>         return super().__str__() + f'\nCourses:\n{self.courses}'
>>  
>>     def assign_course(self, new_course):
>>         self.courses.append(new_course)
>>  
>>     def unassign_course(self, course):
>>         self.courses.remove(course)
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

### Composition and Aggregation
In addition to primitive types, variables in a class can be instances of different classes.
This can often be useful to group relevant data together within a class, such as the molecule class in the Encapsulation lesson, or because a class needs to have ownership of other objects.
There are two different forms that this can take, Composition and Aggregation.
The main difference is the ownership of the object.

Consider the university example we have been using. A university has a large number of students and faculty, but they are not owned by the university. If the university closes, the students and faculty still exist, they just attend or work for a different university. A university is an aggregation of students and faculty.
A university owns courses, if the university closes then the courses cease to exist. A university is composed of courses.

Currently, the Faculty and Students use courses as strings, let us consider extending the courses into its own class.
~~~
class Course:
    def __init__(self, name, description, prerequisitess):
        self.name = name
        self.description = description
        self.prerequisitess = prerequisitess
~~~
{: .language-python}

Here we have a `Course` that consists of a `name` for the course, a `description` for the course, and a list of `prerequisitess` courses that need to be taken.

Now when we look at our `Faculty` class, we see that each Faculty member has a number of courses assigned to them and we utilize the `Course` class to define what each course is. The Faculty member is using the courses, but does not have ownership of them. If the Faculty member leaves the university, the courses will persist and be passed to a new Faculty member to cover. This is an example of aggregation.
