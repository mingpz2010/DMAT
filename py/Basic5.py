__author__ = "mingpingzhou"

from math import sqrt as foobar
# Maybe illustrate some basic control sentences in Python, such as condition, loop etc.
# 1. Print and Import information
# in Python, ':' is used to identify the beginning of sentences, indent technique follows
def basic_info():
    print 'Age : ', 42
    print foobar(4.0)
    x, y, z = 1, 2, 3
    print x, y, z
    x, y = y, x
    print x, y, z
    values = 1,2,3
    print values  # Tuple ADT
    x, y, z = values
    print x
    scoundrel = {'name':'Robin', 'girlfriend': 'Marion'}
    key, value = scoundrel.popitem()
    print key, value

# 2. Condition sentences
def condition_info():
    print True + False + 42
    print bool("")
    print bool(0)
    name = "ming ping zhou"
    if name.endswith('zhou'):
        print "Hello, Ping Zhou!"
    else:
        print "Hello, Stranger!"
    num = 0
    if num > 0:
        print "The number is positive"
    elif num < 0:
        print "The number is negative"
    else:
        print "The number is zero"

    if 'ming' in name:
        print "Maybe the name is MingPingzhou"
    age = 10
    assert 0 < age < 100, 'The age must be realistic'

# 3. Loop sentences
def loop_info():
    x = 1
    while x <= 100:
        print x
        x += 1
    words = ['this', 'is', 'an', 'ex', 'parrot']
    for word in words:
        print word
    numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    for number in numbers:
        print number
    print range(0, 10)
    for number in range(1, 101):
        print number
    d = {'x': 1, 'y' : 2, 'z': 3}
    for key in d:
        print key, 'corresponds to', d[key]
    for key, value in d.items():
        print key, 'corresponds to', value

def useful_func_info():
    print sorted([4, 3, 6, 8, 3])
    print sorted('Hello, world!')
    print list(reversed('Hello, world!'))
    print ''.join(reversed('Hello, world!'))
    for n in range(99, 0, -1):
        root = foobar(n)
        if root == int(root):
            print n
            break

# List comprehension in Python, this is an important technique.
def list_compreh_demo():
    print [x*x for x in range(10)]
    print [x*x for x in range(10) if x % 3 == 0]
    print [(x,y) for x in range(3) for y in range(3)]
    pass
    print x
    del x

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 5"
    basic_info()
    condition_info()
    loop_info()
    useful_func_info()
    list_compreh_demo()

