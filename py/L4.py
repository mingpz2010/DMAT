__author__ = "mingpingzhou"

import math
# Define the basic function to get the value on the discrete points
def my_func(x):
    return math.sqrt(x)

# Numerical integration and numerical differentiation
# 1. Numerical integration
#
# Average partition of integral interval
def Definite_integral_solver(a, b, n):
    if (n <= 0):
        print "Definite integration solver input parameters error!"
        return -1
    if (a == b):
        return 0

    step = b-a/n
    # Trapezoidal formula
    if (n == 1):
        temp = 0.5 * (my_func(a) + my_func(b)) * (b-a)
        return temp
    # Simpson formula
    if (n == 2):
        temp = ((b-a)/6) * (my_func(a) + my_func(b) + 4 * my_func((a+b)/2))

    # the Cotes formula
    if (n == 4):
        x0 = a
        x1 = a + step
        x2 = a + step * 2
        x3 = a + step * 3
        x4 = a + step * 4
        temp = 7*my_func(x0)+32*my_func(x1)+12*my_func(x2)+32*my_func(x3)+7*my_func(x4)
        temp = ((b-a) * temp) / 90.0
        return temp

# 2. Numerical difference



if __name__ == "__main__":
    print "Here are some typical examples"
    # Compute sqrt(x)| (1,9), n = 4
    result = Definite_integral_solver(1, 9, 4)
    print "THE INTEGRATION OF sqrt(x) IS %.10le !" % (result)
