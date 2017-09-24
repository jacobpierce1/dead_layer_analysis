# this script is used to test the speed of class function allocator vs regular allocator.


# test speed of a second classmethod 'constructor'
class test_classmethod_efficiency( object ):
    def __init__( self, x, y ):
        self.x = x
        self.y = y

    @classmethod
    def init2( cls, x, y ):
        return cls( x, y )


# results:
# In [9]: %timeit test_classmethod_efficiency( 1, 2 )
# The slowest run took 10.02 times longer than the fastest. This could mean that an intermediate result is being cached.
# 1000000 loops, best of 3: 405 ns per loop

# In [10]: %timeit test_classmethod_efficiency.init2( 1, 2 )
# The slowest run took 11.48 times longer than the fastest. This could mean that an intermediate result is being cached.
# 1000000 loops, best of 3: 519 ns per loop



# test speed of using a kwarg to specify different construction
class test_kwarg_efficiency( object ):

    def __init__( self, x, y, option=0 ):
        self.x = x
        self.y = y
  

#results:
# In [14]: %timeit test_kwarg_efficiency( 1, 2, option=0 )
# The slowest run took 21.52 times longer than the fastest. This could mean that an intermediate result is being cached.
# 1000000 loops, best of 3: 510 ns per loop

# In [15]: 

# In [15]: 

# In [15]: %timeit test_kwarg_efficiency( 1, 2, option=1 )
# The slowest run took 13.44 times longer than the fastest. This could mean that an intermediate result is being cached.
# 1000000 loops, best of 3: 532 ns per loop
    



# test speed of using a subclass as a different constructor

class test_subclass_constructor_efficiency( test_kwarg_efficiency ):

    def __init__( self, x, y ):
        self.x = x
        self.y = y

# results:
# In [18]: %timeit test_subclass_constructor_efficiency( 1, 2 )
# The slowest run took 17.49 times longer than the fastest. This could mean that an intermediate result is being cached.
# 1000000 loops, best of 3: 409 ns per loop



# conclusion: you can give a class a __init__ constructor that has a
# kwarg or typecheck and can increase efficiency to nearly the same
# efficiency as not having either a kwarg or a typecheck by defining a
# subclass with nothing different except a barebones constructor and
# using that for trustworthy computations within the class. this is
# the approach i have decided to take in meas.py
