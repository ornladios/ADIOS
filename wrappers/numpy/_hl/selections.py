import numpy as np

def select(shape, args):
    """
    Return either SimpleSelection or PointSelection
    """
    if not isinstance(args, tuple):
        args = (args,)

    for a in args:
        if isinstance(a, slice) and a.step is not None and a.step != 1:
            sel = PointSelection(shape)
            sel[args]
            return sel
        if not isinstance(a, slice) and a is not Ellipsis:
            try:
                int(a)
            except Exception:
                sel = PointSelection(shape)
                sel[args]
                return sel

    sel = SimpleSelection(shape)
    sel[args]
    return sel

class Selection(object):
    """
    Abstact selection
    """

    @property
    def shape(self):
        """ Shape of whole dataspace """
        return self._shape

    @property
    def mshape(self):
        """ Shape of current selection """
        raise NotImplementedError()

    def __init__(self, shape):
        shape = tuple(shape)
        self._shape = shape

class SimpleSelection(Selection):
    """
    Simple boundingbox selection
    """

    @property
    def mshape(self):
        """ Shape of current selection """
        return self._mshape

    @property
    def sel(self):
        """ Shape of current selection """
        return self._sel

    def __init__(self, shape, *args, **kwds):
        Selection.__init__(self, shape, *args, **kwds)
        rank = len(self.shape)
        self._sel = ((0,)*rank, self.shape, (1,)*rank, (False,)*rank)
        self._mshape = self.shape

    def __getitem__(self, args):

        if not isinstance(args, tuple):
            args = (args,)

        if self.shape == ():
            if len(args) > 0 and args[0] not in (Ellipsis, ()):
                raise TypeError("Invalid index for scalar dataset (only ..., () allowed)")
            return self

        start, count, step, scalar = _handle_simple(self.shape,args)

        self._sel = (start, count, step, scalar)

        self._mshape = tuple(x for x, y in zip(count, scalar) if not y)

        return self

class PointSelection(Selection):
    """
    Point selection
    """

    @property
    def mshape(self):
        """ Shape of current selection """
        return self._mshape

    @property
    def sel(self):
        """ Shape of current selection """
        return self._sel

    def __init__(self, shape, *args, **kwds):
        Selection.__init__(self, shape, *args, **kwds)
        rank = len(self.shape)

    def __getitem__(self, args):
        raise NotImplementedError("Not implemented yet")

def _expand_ellipsis(args, rank):
    """ Expand ellipsis objects and fill in missing axes.
    """
    n_el = sum(1 for arg in args if arg is Ellipsis)
    if n_el > 1:
        raise ValueError("Only one ellipsis may be used.")
    elif n_el == 0 and len(args) != rank:
        args = args + (Ellipsis,)

    final_args = []
    n_args = len(args)
    for idx, arg in enumerate(args):

        if arg is Ellipsis:
            final_args.extend( (slice(None,None,None),)*(rank-n_args+1) )
        else:
            final_args.append(arg)

    if len(final_args) > rank:
        raise TypeError("Argument sequence too long")

    return final_args

def _handle_simple(shape, args):
    """ Process a "simple" selection tuple, containing only slices and
        integer objects.  Return is a 4-tuple with tuples for start,
        count, step, and a flag which tells if the axis is a "scalar"
        selection (indexed by an integer).

        If "args" is shorter than "shape", the remaining axes are fully
        selected.
    """
    args = _expand_ellipsis(args, len(shape))

    start = []
    count = []
    step  = []
    scalar = []

    for arg, length in zip(args, shape):
        if isinstance(arg, slice):
            x,y,z = _translate_slice(arg, length)
            s = False
        else:
            try:
                x,y,z = _translate_int(int(arg), length)
                s = True
            except TypeError:
                raise TypeError('Illegal index "%s" (must be a slice or number)' % arg)
        start.append(x)
        count.append(y)
        step.append(z)
        scalar.append(s)

    return tuple(start), tuple(count), tuple(step), tuple(scalar)

def _translate_int(exp, length):
    """ Given an integer index, return a 3-tuple
        (start, count, step)
        for hyperslab selection
    """
    if exp < 0:
        exp = length+exp

    if not 0<=exp<length:
        raise ValueError("Index (%s) out of range (0-%s)" % (exp, length-1))

    return exp, 1, 1

def _translate_slice(exp, length):
    """ Given a slice object, return a 3-tuple
        (start, count, step)
        for use with the hyperslab selection routines
    """
    start, stop, step = exp.indices(length)
        # Now if step > 0, then start and stop are in [0, length];
        # if step < 0, they are in [-1, length - 1] (Python 2.6b2 and later;
        # Python issue 3004).

    if step < 1:
        raise ValueError("Step must be >= 1 (got %d)" % step)
    if stop < start:
        raise ValueError("Reverse-order selections are not allowed")

    count = 1 + (stop - start - 1) // step

    return start, count, step
