__author__ = 'mingpingzhou'

# Abstract: basic concepts of computing mathematics.

def error_explain():
    a = 0.6321
    print ("n = %d, a = %.8lf"%(0, a))
    for i in range(1, 10):
        a = 1 - i * a
        print ("n = %d, a = %.8lf"%(i, a))

def valid_explain():
    a = 0.0684
    print ("n = %d, a = %.8lf"%(9, a))
    for i in range(8, -1, -1):
        a = (1 - a) / ( i + 1)
        print ("n = %d, a = %.8lf"%(i, a))

# P1. Avoid x/y, while |y|<<|x|
def Principle_1(x, y, epsilon):
    if abs(y) == 0:
        print "P1 DIV ZERO ERROR!"
        return 0
    else:
        result = abs(x)/abs(y)
        if result<=epsilon:
            print "P1 ERROR!"
            return 0
        return result

# P2. Avoid x - y, the number of x is quite near to y
def Principle_2(x, y, epsilon):
    result = x-y
    if abs(result) <= epsilon:
        print "P2 ERROR!"
        return 0
    return result

# P3. Pay attention to the situation of bigger number vs. small numnber

# P4. Decrease arithmetic times through algorithms' analysis

if __name__ == "__main__":
    print "start to expalin ERROR"
    error_explain()
    print "start to explain VALID"
    valid_explain()

