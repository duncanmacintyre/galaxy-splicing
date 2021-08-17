import numpy as np

# like numpy.where() but can take multiple conditions
#
#       whereswitch(cond1, x1, [cond2, x2, [cond3, x3, ...]][, y], force_object=False)
#
# returns array with values from x1 where cond1, x2 where cond2 but not cond1,
# x3 where cond3 but neither cond1 or cond2, and so on
#
# values from y are chosen where none of cond1, cond2, ... are true
# y defaults to numpy.nan if not given
#
# set force_object=True to ensure that the output has dtype object; this is useful when dealing with
# strings to make sure that numpy.nan does not get cast as a string
#
# whereswitch(cond1, x1, [cond2, x2, ...], y)
#   is equivalent to
# np.where(cond1, x1, np.where(cond2, x2, np.where(..., y)))
#
# e.g.
#   >>> import numpy as np
#   >>> x = np.array([1, 2, 3, 4, 5])
#   >>> y = np.array([6, 7, 8, 9, 10])
#   >>> whereswitch(x < 3, y)
#           numpy.array([6, 7, NaN, NaN, NaN])
#   >>> whereswitch(x < 3, y, x)
#           numpy.array([6, 7, 3, 4, 5])
#   >>> whereswitch(x < 3, x, x < 5, y)
#           numpy.array([1, 2, 8, 9, NaN])
#   >>> whereswitch(x < 3, x, x < 5, y, x < 5, x, 'hello')
#           numpy.array([1, 2, 8, 9, 'hello'])
#   >>> whereswitch(x < 3, 'hello', x < 5, y, x < 5, x, np.nan)                                                                
#           array(['hello', 'hello', '8.0', '9.0', 'nan'], dtype='<U32')
#   >>> whereswitch(x < 3, 'hello', x < 5, y, x < 5, x, np.nan, force_object=True)                                             
#           array(['hello', 'hello', 8, 9, nan], dtype=object)
#

def whereswitch(*args, force_object=False):

    assert len(args) > 1 # we require at least two arguments
    # if force_object, convert all even-numbered inputs to type object first
    if force_object:
        args = [np.asarray(a, dtype='object') if n%2 != 0 else a for n, a in enumerate(args)]
    return _whereswitch(args)

# helper function for recursion
def _whereswitch(args):
    if len(args) == 2: # base case: two arguments - call np.where with np.nan as value to fill in gaps
        return(np.where(*args, np.nan))

    elif len(args) == 3: # base case: three arguments - mimic np.where behaviour
        return(np.where(*args))

    else: # recursion case: four or more arguments
        return np.where(args[0], args[1], _whereswitch(args[2:]))

#   >>> a, b, c = multiwhereswitch(x < 3, ['hello', 'goodbye', np.nan],
#                                  x < 5, [y, x, 'wow'],
#                                  force_object=[True, False, True])
#   >>> a                                                                                                                     
#           array(['hello', 'hello', 8, 9, nan], dtype=object)
#   >>> b                                                                                                                     
#           array(['goodbye', 'goodbye', '3.0', '4.0', 'nan'], dtype='<U32')
#   >>> c                                                                                                                     
#           array([nan, nan, 'wow', 'wow', nan], dtype=object)

def multiwhereswitch(*args, force_object=False):

    return(*[whereswitch(*[a[k] if m%2 != 0 else a for m, a in enumerate(args)],
                         force_object=force_object if type(force_object) == type(True) else force_object[k])
           for k in range(len(args[1]))],)


