__author__ = "mingpingzhou"

import exceptions

# Exception process:
# Exception: base class of all exceptions
# AttributeError: = failure
# IOError: no file open
# IndexError: index is not exist
# KeyError: has no mapping key
# NameError: no names or variables are found
# SyntaxError: codes are error style
# TypeError: wrong type variables
# ValueError: wrong fit values
# ZeroDivisionError: div 0 error!
class MuffledCalculator:
    muffled = False
    def calc(self, expr):
        try:
            return eval(expr)
        except ZeroDivisionError:
            if self.muffled:
                print 'Division by zero is illegal'
            else:
                raise

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 8"
    try:
        x = 10
        y = 0
        print x/y
    except ZeroDivisionError:
        print "The y number can't be zero!"
    calculator = MuffledCalculator()
    # calculator.calc('10/2')
    # calculator.calc('10/0')
    calculator.muffled = True
    calculator.calc('10/0')

