__author__ = "mingpingzhou"

# abstraction, start to write longer program by python.
# sometimes, abstraction is just to implement various functions
def hello(name):
    return 'Hello, ' + name + '!'

def fibs(num):
    result = [0, 1]
    for i in range(num-2):
        result.append(result[-2] + result[-1])
    return result

def square(x):
    'Calculates the square of the number x.'
    return x * x

# The formal parameters, such as string, number and tuple, can't be changed
# But list ,sequence or dictionary can be changed when they are formal parameters
def init(data):
    data['first'] = {}
    data['middle'] = {}
    data['last'] = {}

def add(x, y):
    return x + y

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 6"
    print hello('Ming')
    print fibs(10)
    print square.__doc__
    storage = {}
    init(storage)
    print storage

# some userful functions
# map()
# filter()
# reduce()
# sum()
# apply()

