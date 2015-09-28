__author__ = "mingpingzhou"

from copy import deepcopy

# This part is about dictionary concept in Python
# I think dictionary is the implementation of HASH structure in Computer Science!!!
# If you use mapping, that is dictionary in Python.
# !!! dict function is not a real function , it is a data type as list, tuple or str.
def dictionary_explain():
    # This is a dictionary.
    phonebook = {'Alice': '2341', 'Beth': '9102', 'Cecil': '3258'}
    print phonebook
    items = [('name', 'Gumby'),('age', 42)]
    d = dict(items)
    print d
    print d['name']
    d = dict(name='Gumby', age=42)
    print d

# Sometimes dictionary is quite the same as sequence, but there are still some special characteristics.
def dict_demo():
    phonebook = {'Alice': '2341', 'Beth': '9102', 'Cecil': '3258'}
    x = {}
    x[42] = 'Foobar'
    print x
    print "Cecil's phone number is %(Cecil)s." % phonebook

def dict_demo2():
    template = "'<html>" \
               "<head><title>%(title)s</title></head>" \
               "<body>" \
               "<h1>%(title)s</h1>" \
               "<p>%(text)s</p>" \
               "</body>'"
    print template
    data = {'title': 'My Home Page', 'text': 'Welcome to my home page!'}
    print template % data

def dict_basic_operations():
    phonebook = {'Alice': '2341', 'Beth': '9102', 'Cecil': '3258'}
    print phonebook.clear()
    phonebook = {'Alice': '2341', 'Beth': ['9102','1111'], 'Cecil': '3258'}
    # shallow copy
    phonebook_copy = phonebook.copy()
    print phonebook_copy
    phonebook_copy['Alice'] = '1221'
    print phonebook
    print phonebook_copy
    phonebook_copy['Beth'].remove('1111')
    print phonebook
    print phonebook_copy
    # deep copy
    d = {}
    print d
    d['names'] = ['Alfred', 'Bertrand']
    print d
    c = d.copy()
    print c
    dc = deepcopy(d)
    d['names'].append('Clive')
    print c
    print dc
    print phonebook.items()
    d = {'x': 1, 'y': 2}
    print d.pop('x')
    print d
    print d.values()


if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 4"
    dictionary_explain()
    dict_demo()
    dict_demo2()
    dict_basic_operations()

