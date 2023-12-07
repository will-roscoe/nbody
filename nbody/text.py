from __future__ import annotations
from .base import Iterable, NumType, Vector
from pint import Quantity, UnitRegistry
import os
from . import errors as e
from tqdm import tqdm
import pandas as pd
ur = UnitRegistry()
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
                  'period','pos','vel','acc','time'), vector_pos=True, vector_vel=False, vector_acc = False, engine=None, plotskip=0, c_mass=0):
        
        self.par = {'raw':output_raw, 'ps': plotskip, 'cm': c_mass}
        self.items = items
        self.target = [None, 0]
        self.engine = engine
        self._m = {'pos':vector_pos, 'vel':vector_vel, 'acc':vector_acc}
        self.q = {'identity':'','mass':0*ur.kg,'radius':0*ur.m,'energy':0*ur.joule,
                  'period':0*ur.s,'pos':0*ur.m,'vel':0*ur.m/ur.s,'acc':0*ur.m/ur.s**2}
    def _lookup(self):
        return {'identity' : self.target[0].identity,
                    'mass'   : self.target[0].mass.c()*ur(self.target[0].mass.units),
                    'radius' : self.target[0].radius.c()*ur(self.target[0].radius.units),
                    'energy' : self.target[0].get_('ke', self.target[1], self.par['ps'], self.par['cm']) * ur.joule,
                    'period' : self.target[0].get_('period', self.target[1], self.par['ps'], self.par['cm'], engine=self.engine) * ur.s,
                    'pos' : ([self.target[0].pos[self.target[1]][n]*ur(self.target[0].pos.units)
                            for n in range(3)] if self._m['pos']
                            else Vector(self.target[0].pos[self.target[1]]).magnitude()*ur(self.target[0].pos.units)),
                    'vel' : ([self.target[0].vel[self.target[1]][n]*ur(self.target[0].vel.units)
                            for n in range(3)] if self._m['vel'] 
                            else Vector(self.target[0].vel[self.target[1]]).magnitude()*ur(self.target[0].vel.units)),
                    'acc' : ([self.target[0].acc[self.target[1]][n]*ur(self.target[0].acc.units) 
                            for n in range(3)] if self._m['acc'] 
                            else Vector(self.target[0].acc[self.target[1]]).magnitude()*ur(self.target[0].acc.units)),
                    'time' : self.target[1]*self.engine.dt * ur.s,
    }
   
    def _basetemplate(self):
        ret = ''
        for key,st in {'identity':"{identity}",
                       'mass': "\nMass: {mass}",
                       'radius':"\nRadius: {radius}",
                       'energy':"\nKE: {energy}",
                       'period':"\nPeriod: {period}",
                       'pos':"\nPOS: {pos}",
                       'vel':"\nVEL: {vel}",
                       'acc':"\nACC: {acc}",
                       'time': "\nT: {time}"}:
            if key in self.items:
                ret+=st
    
    def _quantities(self, body=None):
        if body != None: self.target[0] = body
        q = {}
        for key in self.items:
            self.q[key] = self._lookup()[key]
        return self.q
    
   
    def _get_best_unit(self, quant, units):
        if units is None or (isinstance(quant, Quantity) and quant.m=='NaN'):
            return quant
        if units[0] == 'check_mass':
            if self.q['mass'] > 0.001*ur.emass:
                return quant
            else:
                units = units[1]
        for _un in units:
            if quant > 0.05*_un:
                return quant.to(_un)
        return quant
        
    
    def convert(self, arg=None):
        conv = {}
        if not arg:
            args = self._quantities()
        else:
            args = arg
        for key,value in args.items():
            if not isinstance(value, Iterable): 
                conv[key] = self._get_best_unit(value, _units[key])
            else:
                conv[key] = [self._get_best_unit(v, _units[key]) for v in value]
        return conv
    def _quantities_to_strings(self, arg=None):
        strings = []
        if arg == None:
            args = self._quantities()
        else:
            args = arg
        for key, value in args.items():
            try:    
                if isinstance(value, Iterable):
                    strings.append('('+''.join((f'{q:.4f~P}' for q in value))+')') #q.to_compact()
                elif isinstance(value, Quantity) and value.m != 'NaN':
                    strings.append(f'{value:.4f~P}') #value.to_compact()
                elif isinstance(value, Quantity) and value.m == 'NaN':
                    strings.append('NaN')
                else:
                    strings.append(value)
            except ValueError:
                strings.append('NaN')
        return strings
    
    def __str__(self):
        if self.par['raw'] == False:
            return self._basetemplate().format(*self._quantities_to_strings(self.convert()))
        else:
            return self._basetemplate().format(*self._quantities_to_strings())
        
def body_from(object,input_type='init'):
    if input_type == 'init':    
        if isinstance(object, str): 
            if os.path.isfile(object):
                params = dict()
                with open(object, 'r') as f:
                    for line in f:
                        p = list(l.strip() for l in line.strip().split('='))
                        if p[0] in ('name', 'identity', 'id'):
                            params['identity'] = p[1]
                        elif p[0] in ('color', 'colour'):
                            params['color'] = p[1]
                        elif p[0] in ('mass', 'radius', 'bounce'):
                            params[p[0]] = float(p[1])
                        elif 'pos' in p[0]:
                            tp = p[1].strip('()[]').split(',')
                            params['init_pos'] = list(float(t) for t in tp)
                        elif 'vel' in p[0]:
                            tp = p[1].strip('()[]').split(',')
                            params['init_vel'] = list(float(t) for t in tp)
                        else:
                            tqdm.write(f'«body_from()» → Couldn\'t parse "{line}" so it has been skipped.')
            else:
                e.raise_fnf_error(object)
                return Body(**params)
        if isinstance(object, dict):
            return Body(**object)
        else:
            e.raise_type_error('object', (dict, str), object)
    
    elif input_type.startswith() == 'pos':
        if os.path.isfile(object):
            result = {}
            with open(object, 'r') as file:
                lines = file.readlines()
                parsing_position = False
                found_dt = False
                positions = []
                dt = 0
                for line in lines:
                    parts = line.split('=')
                    if len(parts) == 2:
                        key, value = map(str.strip, parts)
                        value = value.replace('(', '').replace(')', '').replace(',', '')
                        if key in ('name', 'identity', 'id'):
                            result['identity'] = value
                        elif key in ('color', 'colour'):
                            result['color'] = value
                        elif key in ('mass', 'radius', 'bounce'):
                            result[key] = float(value)
                        elif key == 'position':
                            parsing_position = True
                            continue
                        elif key == 'dt':
                            found_dt = True
                            dt = float(value)

                    if parsing_position:
                        if line.strip() == ')':
                            result['position'] = positions
                            parsing_position = False
                        else:
                            position_str = line.strip().replace('(', '').replace(')', '')
                            position_tuple = tuple(map(int, position_str.split(',')))
                            positions.append(position_tuple)
            if found_dt == True:
                params = result
                params['init_pos'] = result['position'][0]
                body = Body(**params)
                for i,x in enumerate(result['position'][1:]):
                    body.pos.next(x)
                    body.vel.next((Vector(result['position'][i])-x)/dt)
                    body.acc.next((0,0,0))
                return body
            else:
                raise LookupError('could not find a dt value.')
        else:
            e.raise_fnf_error(object)
    else:
        e.raise_value_error('input_type', str, input_type)

def engine_from(object):   
    if isinstance(object, str):
        if os.path.isfile(object):
            with open(object, 'r') as f:
                _input = f.readlines()
                # getting parameters for engine
                eng_params = dict()
                for line in _input:
                    parts = line.split('=')
                    if len(parts) == 2:
                        key, value = map(str.strip, parts)
                        value = value.replace('(', '').replace(')', '').replace(',', '')
                        if key == 'dt':
                            eng_params['dt'] = float(value)
                        elif key == 'checking_range':
                            eng_params['checking_range'] = float(value)
                # getting start and end indexes of each body
                body_strings = dict()
                func_strings = []
                for i,line in enumerate(_input):
                    if line.startswith('*'):
                        key = line.strip('[\n').strip(']').strip('*').strip()
                        body_strings[key] = i
                        for s_i,s_line in enumerate(_input):
                            if s_i > i:
                                if ']' in s_line and isinstance(body_strings[key], NumType):
                                    body_strings[key] = _input[i:s_i]
                    elif line.startswith('-*'):
                        val = [p.strip() for p in line.strip('-*').strip('\n').split('=')]
                        func_strings.append(val)
                bodies = dict()
                for key,value in body_strings.items():
                    params = dict()
                    for line in value:
                        p = list(l.strip() for l in line.strip().split('='))
                        if p[0] in ('name', 'identity', 'id'):
                            params['identity'] = p[1]
                        elif p[0] in ('color', 'colour'):
                            params['color'] = p[1]
                        elif p[0] in ('mass', 'radius', 'bounce'):
                            params[p[0]] = float(p[1])
                        elif p[0] in ('pos', 'position'):
                            tp = p[1].strip('()[]').split(',')
                            params['init_pos'] = list(float(t) for t in tp)
                        elif p[0] in ('vel', 'velocity'):
                            tp = p[1].strip('()[]').split(',')
                            params['init_vel'] = list(float(t) for t in tp)
                        elif '[' in p[0] or ']' in p[0] or '*' in p[0] or '#' in p[0]:
                            pass
                        else:
                            tqdm.write(f'«engine_from()» → Couldn\'t parse "{line}" on line {i} so it has been skipped.')
                    bodies[key] = Body(**params)
                engine = Engine(**eng_params)
                engine.attach_bodies(list(bodies.values()))
                for val in func_strings:
                    if val[0] == 'plane' and len(val) == 3:
                        engine.create_plane(str(val[1]), float(val[2]))
                    elif val[0] == 'field' and len(parts) == 2:
                        vect = [float(v) for v in val[1].replace('(', '').replace(')', '').split(',')]
                        engine.create_acceleration(vect)
                    elif val[0] == 'rel' and len(val) == 2:
                        engine.make_relative_to(bodies[val[1]])
                    elif val[0] == 'sim' and len(val) == 2:
                            engine.simulate(int(val[1]))
                    else:
                        tqdm.write(f'«engine_from()» → Couldn\'t parse function "{' '.join(val)}" so it has been skipped.')
                return engine
        else:
            e.raise_fnf_error(object)
    else:
        e.raise_type_error('object', str, object)

def export_data(object, loc, overwrite=True):
    if isinstance(object, Engine):
        os.makedirs(loc, exist_ok=overwrite)
        eng_info = [
            f'dt = {object.dt}\n',
            f'checking_range = {object._rangechk}\n',
            f'number of intervals: {len(object)}\n',
            f'# bodies:({len(object.bodies)})\n',
        ]
        ids = []
        for b in object.bodies:  
            obj_id = str(b.identity).replace(' ','_').split(':')[0].split('(Static)')[0]
            ids.append(obj_id)
            for line in [f'*{obj_id} [\n',
                         f'name = {b.identity}\n',
                         f'mass = {b.mass.c()}\n',
                         f'radius = {b.radius.c()}\n',
                         f'color = {b.color}\n',
                         f'bounce = {b.bounce}\n',
                         f'position = {b.pos.c()}\n',
                         f'velocity = {b.vel.c()}\n' 
                         ']\n']:
                eng_info.append(line)
        with open(f'{loc}/eng_info.txt', ('w' if overwrite is True else 'x')) as eng:
            eng.writelines(eng_info)
        for i,b in enumerate(object.bodies):
                info_dict = {'dt':[t*object.dt for t in range(len(object))],
                             'posX':b.pos.X.record,
                             'posY':b.pos.Y.record,
                             'posZ':b.pos.Z.record,
                             'velX':b.vel.X.record,
                             'velY':b.vel.Y.record,
                             'velZ':b.vel.Z.record,
                             'accX':b.acc.X.record,
                             'accY':b.acc.Y.record,
                             'accZ':b.acc.Z.record
                             }
                d_df = pd.DataFrame.from_dict(info_dict, 'index').transpose()
                d_df.to_csv(f'{loc}/{ids[i]}.csv', sep=',',index=False, mode=('w' if overwrite is True else 'x'))





from .core import Body, Engine


