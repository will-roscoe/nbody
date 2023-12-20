
from pint import Quantity, UnitRegistry
# ⤤ unit parsing
from ..core.base import Vector, Iterable, _O
ur = UnitRegistry()
Q_ = ur.Quantity
ur.define('light_year = 9460730472580800 * meter = ly = lightyear')
ur.define('jupiter_mass = 1.89818*10**27 * kg = $M_J$ = jmass')
ur.define('solar_mass = 1.9885*10**30 * kg = $M_O$ = smass')
ur.define('earth_mass = 5.97219*10**24 * kg = $M_E$ = emass')
ur.define('astronomical_unit_per_year = astronomical unit / year')

### conversion preferences:
## (in order of largest to smallest)
mass_conversions = (ur.smass,
                    ur.jmass,
                    ur.emass)

distance_conversions = (ur.light_year,
                        ur.astronomical_unit)

velocity_conversions = (ur.astronomical_unit/ur.day,
                        ur.km/ur.second,
                        ur.km/ur.hour,
                        ur.km/ur.day)

time_conversions = (ur.year,ur.month,ur.day)

acceleration_conversions = (ur.km/(ur.second**2),
                            ur.km/(ur.hour**2),
                            ur.km/(ur.day**2))

####

_units = {'identity': None,
        'mass':mass_conversions,
        'radius':distance_conversions,
        'energy': None,
        'period':time_conversions,
        'pos':distance_conversions,
        'vel':('check_mass', velocity_conversions),
        'acc':('check_mass', acceleration_conversions),
        'time':time_conversions}
ur.default_format = '~P'

class Formatter:
    def __init__(self, output_raw = False, items=('identity','mass','radius','energy',
                  'period','pos','vel','acc','time'), 
                  vector_pos=True, vector_vel=False, vector_acc = False, engine=None, plotskip=0, c_mass=0):
        # sort parameters into useful form
        self.par = {'raw':output_raw, 'ps': plotskip, 'cm': c_mass}
        self.items = items
        self.target = [None, 0] # target body and index
        self.engine = engine
        # whether to do mag or vector
        self._m = {'pos':vector_pos, 'vel':vector_vel, 'acc':vector_acc}
        # base quantity dict
        self.q = {'identity':'','mass':0*ur.kg,'radius':0*ur.m,'energy':0*ur.joule,
                  'period':0*ur.s,'pos':0*ur.m,'vel':0*ur.m/ur.s,'acc':0*ur.m/ur.s**2}
    
    def _lookup(self):
        # dict containing pint quantities of target
        return {
    'identity' : self.target[0].identity,
    'mass'   : Q_(float(self.target[0].mass.c()), ur(self.target[0].mass.units)),
    'radius' : Q_(float(self.target[0].radius.c()), ur(self.target[0].radius.units)),
    'energy' :Q_(float(_O(self.target[0].get_('ke', self.target[1], self.par['ps'], self.par['cm']))), ur.joule),
    'period' : Q_(float(_O(self.target[0].get_('period', self.target[1], self.par['ps'], self.par['cm'], 
                                               engine=self.engine))), ur.s),
    'pos' : ([Q_(float(self.target[0].pos[self.target[1]][n]),ur(self.target[0].pos.units))
            for n in range(3)] if self._m['pos']
            else Q_(float(_O(Vector(self.target[0].pos[self.target[1]]).magnitude())), ur(self.target[0].pos.units))),
    'vel' : ([Q_(float(self.target[0].vel[self.target[1]][n]), ur(self.target[0].vel.units))
            for n in range(3)] if self._m['vel'] 
            else Q_(_O(float(Vector(self.target[0].vel[self.target[1]]).magnitude())), ur(self.target[0].vel.units))),
    'acc' : ([Q_(float(self.target[0].acc[self.target[1]][n]), ur(self.target[0].acc.units)) 
            for n in range(3)] if self._m['acc'] 
            else Q_(float(_O(Vector(self.target[0].acc[self.target[1]]).magnitude())), ur(self.target[0].acc.units))),
    'time' : Q_(self.target[1]*self.engine.dt, ur.s),
    }
    
    def _basetemplate(self):
        # template string with formattable entries
        return ''.join([{'identity':"{identity}",
                       'mass': "\nMass: {mass}",
                       'radius':"\nRadius: {radius}",
                       'energy':"\nKE: {energy}",
                       'period':"\nPeriod: {period}",
                       'pos':"\nPOS: {pos}",
                       'vel':"\nVEL: {vel}",
                       'acc':"\nACC: {acc}",
                       'time': "\nT: {time}"}[i] for i in self.items])
    
    def _quantities(self, body=None):
        if body != None: 
            self.target[0] = body
        q = {}
        for key in self.items:
            q[key] = self._lookup()[key]
        return q
    
    
    def _get_best_unit(self, quant, units):
        if units is None or (isinstance(quant, Quantity) and quant.m=='NaN'):
            return quant
        # if we want to leave as is if object is small enough
        if units[0] == 'check_mass':
            if self.q['mass'] > 0.001*ur.emass:
                return quant
            else:
                units = units[1]
        # test for best unit.
        for _un in units:
            if quant > 0.05*_un:
                return quant.to(_un)
        return quant
        
    
    def convert(self, arg=None):
        conv = {}
        args = arg if arg else self._quantities()
        # get best units for each item we want
        for key,value in args.items():
            if not isinstance(value, (list, tuple)): 
                # scalars
                conv[key] = self._get_best_unit(value, _units[key])
            else:
                # vectors : find best for each.
                conv[key] = [self._get_best_unit(v, _units[key]) for v in value]
        return conv
    
    def _quantities_to_strings(self, arg=None):
        strings = {}
        args = self._quantities() if arg == None else arg
        for key, value in args.items():
            try:    
                if isinstance(value, Iterable):
                    # vectors
                    strings[key] = ('('+''.join((f'{q:.4f~P}' for q in value))+')') #q.to_compact()
                elif isinstance(value, Quantity) and value.m != 'NaN':
                    # scalars
                    strings[key] = (f'{value:.4f~P}') #value.to_compact()
                elif isinstance(value, Quantity) and value.m == 'NaN':
                    # NaN values (not neccesarily a str)
                    strings[key] = ('NaN')
                else:
                    # string values
                    strings[key] = (value)
            except ValueError:
                strings[key] = ('NaN')
        return strings
    
    def __str__(self):
        if self.target[0] is not None:
            if self.par['raw'] == False:
                # if called as a string, return _basetemplate with the correct items inserted.
                return self._basetemplate().format(**self._quantities_to_strings(self.convert()))
            else:
                return self._basetemplate().format(*self._quantities_to_strings())
        else:
            return ''

globals