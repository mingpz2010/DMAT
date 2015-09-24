__author__ = "mingpingzhou"

# Learn how to store data in complex structure such as list and dictionary
# I think this part is the complex python program's basic, so it should be mastered.

# Sequence ADT
# Basic operations: indexing, sliceing, adding, multiplying
# When sliceing, the first element is contained but the last isn't
def sequence_explain():
    edward = ['Edward Gumby', 42]
    john = ['John Smith', 50]
    database = [edward, john]
    print database
    greeting = "Hello"
    print greeting[0]
    print greeting[-1]
    months = ['January','February','March','April','May','June',
              'July','August','September','October','November','December']
    endings = ['st','nd','rd']+17*['th']+['st','nd','rd']+7*['th']+['st']
    year = "1985"
    month = "9"
    day = "21"
    month_name = months[int(month)-1]
    ordinal = day + endings[int(day)-1]
    print month_name + ' ' + ordinal + ', ' + year
    tag = '<a href="http://www.python.org">Python web site</a>'
    print tag[9:30]
    print tag[32:-4]

def url_split(url) :
    domain = url[11:-4]
    print "Domain name : " + domain

def sequence_arithmetic():
    a = [1,2,3]
    b = [4,5,6]
    print a + b
    print a * 5
    print [None] * 5
    sentence = "Python is very interesting!"
    screen_width = 80
    text_width = len(sentence)
    box_width = text_width + 6
    left_margin = (screen_width - box_width) / 2
    print ' '
    print ' '*left_margin+'+'+'-'*(box_width-2)+'+'
    print ' '*left_margin+'| '+' '*text_width+' |'
    print ' '*left_margin+'| '+ sentence +' |'
    print ' '*left_margin+'| '+' '*text_width+' |'
    print ' '*left_margin+'+'+'-'*(box_width-2)+'+'
    print ' '
    permissions = 'rw'
    print 'w' in permissions
    print 'x' in permissions
    numbers = [100,34,678]
    print len(numbers)
    print max(numbers)
    print min(numbers)

# List ADT, sometimes the list is mutable
# Basic operations: the same as sequence ADT, but it has some operations that can modify list, such as delete etc.
def list_explain():
    print list('Hello')
    x = [1, 1, 1]
    x[1] = 2
    print x
    names = ['Alice','Beth','Cecil','Dee-Dee','Earl']
    del names[2]
    print names
    name = list('Perl')
    print name
    name[2:] = list('ar')
    print name

def list_method_demo():
    lst = [1,2,3]
    lst.append(4)
    print lst
    print ['to','be','or','not','to','be'].count('to')
    y = [[1,2], 1, 3, [2,1,[1,2]]]
    print y.count(1)
    print y
    print y.count([1,2])
    print y.index(3)
    print lst.pop()
    print lst
    print lst.reverse()
    print lst
    x = [4,6,2,1,7,9]
    y = x[:]
    y.sort()
    print x
    print y

# Tuple ADT
# Somethimes the tuple ADT is viewed as a whole variable, so it can be recognized as structure.
# So in scientific research, tuple is a data type concept, not a function concept.
def tuple_explain():
    t = (1,2,3,4,5,6)
    print t
    single_t = (42,)
    print single_t
    print tuple([1,2,3,4,5,6])

if __name__ == "__main__":
    print "RUN <Python Basic Tutorial> Contention 2"
    sequence_explain()
    url_split("http://www.something.com")
    sequence_arithmetic()
    list_explain()
    list_method_demo()
    tuple_explain()

# Here list some useful functions:
# cmp()
# len()
# list()
# max()
# min()
# reverse()
# sorted()
# tuple()