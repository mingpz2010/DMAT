__author__ = 'mingpingzhou'

def my_func(x):
    ans = x * x * x - x - 1
    return ans

# solve non-linear equation

# 1. Bisection method is the basic of non-linear equation method!
def bisection_solve(a, b):
    x = (a+b)/2
    x_last = x + 1
    print "X = %.10le" %(x)
    while (abs(x-x_last) >= 1e-6):
        x_last = x
        fans = my_func(x)
        if (fans > 0):
            b = x
        elif (fans < 0):
            a = x
        else:
            break
        x = (a+b)/2
        print "X = %.10le" %(x)
    print "FINAL ANS IS X = %.10le" %(x)

if __name__ == "__main__":
    bisection_solve(1.0, 1.5)

