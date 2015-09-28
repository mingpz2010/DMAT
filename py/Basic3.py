__author__ = "mingpingzhou"

from math import pi
from string import Template

# This module is related to the operation of string in Python
# Remember, the string is a constant variable!!!
def string_expalin():
    format = "Pi with three decimals : %.3f"
    print format % pi
    s = Template('$x. glorious $x!')
    print s.substitute(x='slurm')
    # if the substitute content is part of the noun word, the parameter must be contained by a bracket
    s = Template("It's ${x}tastic!")
    print s.substitute(x='slurm')

def string_output_demo():
    print '%s plus %s equals %s' %(1, 1, 2)
    print '%010.2f' % pi
    print '%-10.2f' % pi
    print ('% 5d' % 10) + '\n' + ('% 5d' % -10)
    print ('%+5d' % 10) + '\n' + ('%+5d' % -10)

def print_priceTable():
    width = 35
    price_width = 10
    item_width = width - price_width
    header_format = '%-*s%*s'
    format = '%-*s%*.2f'
    print '=' * width
    print header_format % (item_width, 'Item', price_width, 'Price')
    print '-' * width
    print format % (item_width, 'Apples', price_width, 0.4)
    print format % (item_width, 'Pears', price_width, 0.5)
    print format % (item_width, 'Cantaloupes', price_width, 1.92)
    print format % (item_width, 'Dried Apricots (16 oz.)', price_width, 8)
    print format % (item_width, 'Prunes (4 lbs.)', price_width, 12)

def other_methods():
    s = 'With a moo-moo here, and a moo-moo there'
    print s.find('moo')
    seq = ['1','2','3','4','5']
    sep = '+'
    print sep.join(seq)
    dirs = '','usr','bin','env'
    print '/'.join(dirs)
    print 'C:' + '\\'.join(dirs)    # Pay attention to reverse character!
    print '/usr/bin/env'.split('/')
    print '         internal whitespace is kept           '.strip()

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 3"
    string_expalin()
    string_output_demo()
    print_priceTable()
    other_methods()

