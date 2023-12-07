from math import prod
import pytest

import nbody.base as nb


ns = [10,0.01,-3.2e8]
vs = [nb.Variable(n) for n in ns]
hvs = [nb.HistoricVariable(n) for n in ns]
def add_sc():
        tval = vs[0]
        tval =+ vs[1] + ns[2]
        assert tval == sum(ns)
def div_sc():
        tval = vs[0]
        tval /= vs[1] / ns[2]
        assert tval == ns[0]/(ns[1]/ns[2])
def mult_sc():
        tval = vs[0]
        tval *= vs[1] * ns[2]
        assert tval == prod(ns)



