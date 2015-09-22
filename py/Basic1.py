__author__ = "mingpingzhou"

import math
import cmath
# or. from math import floor

print 10 + pow(2, 3*5)/3.0
print math.floor(32.9)
print cmath.sqrt(-1)

print repr("Hello World!")
print repr(100000L)
print str("Hello World!")
print str(100000L)

# When use unicode string, please add character 'u' on the head of the string
print u"This is Unicode String!"

def user_input():
    x = input("The meaning of life : ")
    print "LIFE Meaning is %d" %(x)
    if x > 60 : print "Your age is quite old!"

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 1"
    # user_input()

# Here list some useful functions:
# abs()
# cmath.sqrt()
# float()
# help()
# input()
# int()
# long()
# math.ceil()
# math.floor()
# math.sqrt()
# pow()
# raw_input()
# repr()
# round()
# str()

