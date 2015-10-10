__author__ = "mingpingzhou"

# Object-oriented: Polymorphism, Encapsulation, Inheritance
__metaclass__ = type

# self is the difference between methods and functions !
# methods have the self parameter
class Person:
    def setName(self, name):
        self.name = name
    def getName(self):
        return self.name
    def greet(self):
        print "Hello, world! I'm %s." % self.name

class Bird:
    song = 'Squaawk!'
    def sing(self):
        print self.song

# use __ trick to ensure your method is private
# or use _ trick, so import sentence won't import your class method into various modules
class Secretive:
    def __inaccessible(self):
        print "Bet you can't see me..."
    def accessible(self):
        print "The secret message is:"
        self.__inaccessible()

# Don't have a construct function
# Class has local variables access property
class MemberCounter:
    members = 0
    def init(self):
        MemberCounter.members += 1

# super class(base class) definition
class Filter:
    def init(self):
        self.blocked = []
    def filter(self, sequence):
        return [x for x in sequence if x not in self.blocked]

# SPAMFilter is the child class of Filter
class SPAMFilter(Filter):
    def init(self):
        self.blocked = ['SPAM']

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 7"
    foo = Person()
    foo.setName('Mingpingzhou')
    foo.greet()
    f = Filter()
    f.init()
    print f.filter([1,2,3])
    s = SPAMFilter()
    s.init()
    print s.filter(['SPAM', 'SPAM', 'SPAM', 'SPAM', 'eggs', 'bacon', 'SPAM'])
    # Check one class is whether the child class of another class!!!
    print issubclass(SPAMFilter, Filter)
    print issubclass(Filter, SPAMFilter)
    # Check one object is whether the instance of one class
    print isinstance(s, SPAMFilter)
    print isinstance(s, Filter)
    print isinstance(s, str)
    print s.__class__

# write down your description of the program, then choose the noun, verb and adj.
# class ----> noun
# method ----> verb
# property -----> adj

