
import pytest
import mpmath as mp
from ..core import base as nb


# building test objects
_hvar = nb.HistoricVariable(0)
_hvar.next([i for i in range(1,10)])
_hvect = nb.HistoricVector(0,0,0)
for i in range(1,10):
        _hvect.next((i,i,i))

class TestBaseFuncs:
        '''Testing class independent functions
        '''
        def test_ntype(self):
                assert nb._ntype(1,1.1) == mp.mpf and nb._ntype('string') == str
        def test_typecheck(self):
                with pytest.raises(TypeError):
                        nb.typecheck(((123, str)))
                assert nb.typecheck(((2, int))) == 2
        def test_o(self):
                assert nb._O(nb.Variable(1)) == 1 and nb._O(nb.NullVector()) == (0,0,0)
        def test_v(self):
                with pytest.raises((IndexError, TypeError)):
                        nb._V([0])
                assert nb._V(0) == 0 and nb._V((0,0,0)).c() == nb.Vector((0,0,0)).c()


class TestVarVectorBasics:
        '''Testing __init__, length and type definitions
        '''
        def test_variable(self):
                _var = nb.Variable(0)
        def test_histvariable(self):
                _hvar = nb.HistoricVariable(0)
        def test_vector(self):
                _vect = nb.Vector((0,0,0))
                _vect2 = nb.Vector(x=0, y=0, z=0)
        def test_histvector(self):
                _hvect = nb.HistoricVector(0,0,0)
                _hvect2 = nb.HistoricVector(li=(0,0,0))
        def test_len(self):
                var_l = len(nb.Variable(0))
                hvar_l = len(_hvar)
                vect_l = len(nb.Vector((0,0,0)))
                hvect_l = len(_hvect)
                assert var_l == vect_l and hvar_l == hvect_l
        def test_types(self):
                assert (isinstance(nb.Variable(0), nb.VarType) and
                        isinstance(nb.HistoricVariable(0), nb.VarType) and
                        isinstance(nb.Vector((0,0,0)), nb.VectorType) and
                        isinstance(nb.HistoricVector(0,0,0), nb.VectorType) and
                        isinstance(nb.NullVector(), nb.VectorType))

class TestVarMaths:
        '''Testing Basic Math Operators beween 1d variables and python numbers
        '''
        def test_mult(self):
                val = nb.Variable(1) * nb.Variable(2) * 3.
                assert val == 6.
        def test_add(self):
                assert nb.Variable(1.) + 3.1 + nb.Variable(2e0) == 6.1 
        def test_sub(self):
                ans = nb.Variable(1) - 3.1 - nb.Variable(-2)
                assert ans == -0.1
        def test_truediv(self):
                assert nb.Variable(1.1) / nb.Variable(2.2) / 3. == 1/6


class TestVectMaths:
        ''' Testing Math operators between 3-vectors and python tuples
        '''
        def test_add(self):
                test = (nb.Vector([1,0.1,2e0]) + (-1,1,1)) + nb.Vector([0,-5.1,3])
                assert test == (0., -4.,6.)
        def test_sub(self):
                test = (nb.Vector([1,0.1,2e0]) - (-1,1,1)) - nb.Vector([0,-5.1,3])
                assert test == (2.,4.2,-2.)  
        def test_dot(self):
                val1 = nb.Vector([1,0.1,2e0]) * (-1,1,1) 
                val2 = nb.Vector([0,-5.1,3]) * nb.Vector((1,1,1))
                val3 = (2,1.,0.6) * nb.Vector((3,-1,0.2))
                assert val1 == 1.1 and val2 == -2.1 and val3 == 5.12
        def test_mag(self):
                v1 = nb.Vector([1,1,1]).magnitude()
                zero = nb.NullVector().magnitude()
                mult = nb.Vector([4, 13, 16]).magnitude()
                assert v1 == mp.sqrt(3) and zero == 0. and mult == 21.
        def test_unit(self):
                vec = nb.Vector([1.23,4.2, 10000])
                assert vec.unit()*vec.magnitude() == vec   

class TestInteroperability:
        ''' Testing translation between classes and math.
        ''' 
        def test_var_to_hvar(self):
                hvar = nb.HistoricVariable(1.)
                hvar.next(nb.Variable(0.))

                assert (hvar.record[0] == 1.
                        and hvar.record[1] == 0. 
                        and nb.HistoricVariable(1.) + nb.Variable(-1.) == 0.)
        def test_full(self):
                n = {'p':nb.HistoricVector(0.,0.,0.), 
                     'v':nb.HistoricVector(2.,-0.5,0.), 
                     'a':nb.HistoricVector(0.01,-0.03,0.5)}
                time = nb.HistoricVariable(0.)
                for i in range(1,10):
                        time.next(i*abs(n['v'].c(1))) 
                        n['a'].next((n['p'].magnitude()/time, n['v'].magnitude()/time, n['a'].magnitude()/time))
                        n['v'].next((n['a']/time + n['v']))
                        n['p'].next((n['v'].unit()/time+n['p']*n['v'].magnitude()))