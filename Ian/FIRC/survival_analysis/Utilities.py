import numpy as np

VERSION = 1.0

#builds a dictionary in python
#directly ported line for line from Marks code
#can probably be changed or removed when improving the code
def buildhash(hashkeys,hashvals):
    
    soughthash = {}
    
    for j in range(len(hashkeys)):
        key = hashkeys[j]
        soughthash[key] = hashvals[j]
        
    return(soughthash)

# rounds a number somehow. Will need to test this to see what it should return and what it actually returns
def round_num(number):
    
    number = int(number + .5*np.sign(number))
    
    return(number)

# In Marks code this function removes whitespace at the beginning and end of an array. But as that isn't a thing in python it will probably be unecessary.
def trim(array):
    
    j = 0
    while j<len(array):
        #i don't think this function is needed in python
        #so the loop is left empty
        
        j = j +1
        
    return(array)