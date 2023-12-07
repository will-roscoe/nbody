from __future__ import annotations
import os
import pandas as pd
from tqdm import tqdm
from . import errors as e
from ..core.base import Vector, NumType
from ..core.body import Body
from ..core.engine import Engine
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



