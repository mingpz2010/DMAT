__author__ = "mingpingzhou"

# This is the meta class definition
__metaclass__ = type

# Magic methods, properties and iterations
# If there is no __metaclass__ declaration, new style and old style are different. New represents new tyle!
# But __metaclass__ is declared, these definitions are the same.
class NewStyle(object):
    print "New Class"
class OldStyle:
    print "Old Class"

# 1. magic of construct method
class FooBar:
    def __init__(self, value=42):
        self.somevar = value

# 2. re-write general methods and special construct methods
class A:
    def hello(self):
        print "Hello, I'm A."
class B(A):
    def hello(self):
        print "Hello, I'm B."

class Bird:
    def __init__(self):
        self.hungry = True
    def eat(self):
        if self.hungry:
            print 'Aaaah...'
            self.hungry = False
        else:
            print "No, thanks!"

class SongBird(Bird):
    def __init__(self):
        # method 1 : This sentence is very important!
        # Bird.__init__(self)
        # method 2 :
        super(SongBird, self).__init__()
        self.sound = 'Squawk!'
    def sing(self):
        print self.sound

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 9"
    f = FooBar()
    print f.somevar
    sb = SongBird()
    sb.sing()
    sb.eat()
    sb.eat()


