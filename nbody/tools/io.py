from __future__ import annotations
import os
import pandas as pd
from mpmath import mpf
from tqdm import tqdm
from . import errors as e
from ..core.base import Vector, NumType, _O
from ..core.body import Body
from ..core.engine import Engine


def tup(iterable):
    tmp = _O(iterable)
    return str(tuple(str(c) for c in tmp)).replace("'", "")


def obj_from(object, obj='engine'):
    if isinstance(object, str):
        if os.path.isfile(object):
            with open(object, 'r') as f:
                # get all lines
                _input = f.readlines()
                if obj == 'engine':
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
                                eng_params['checking_range'] = int(value)
                    # sorting special items
                    body_strings = dict()
                    func_strings,var_strings = [], []
                    for i,line in enumerate(_input):
                        # body objects: find end of body object definition
                        if line.startswith('*'):
                            key = line.strip('[\n').strip(']').strip('*').strip()
                            body_strings[key] = i
                            for s_i,s_line in enumerate(_input):
                                if s_i > i:
                                    if ']' in s_line and isinstance(body_strings[key], NumType):
                                        body_strings[key] = _input[i:s_i]
                        # callables
                        elif line.startswith('-*'):
                            val = [p.strip() for p in line.strip('-*').strip('\n').split('=')]
                            func_strings.append(val)
                        # attribute definitions
                        elif line.startswith('~'):
                            val = [p.strip() for p in line.strip('~').strip('\n').split('=')]
                            var_strings.append(val)
                    bodies = dict()
                    # get body info and create each body
                    for key,value in body_strings.items():
                        params = dict()
                        for line in value:
                            p = list(l.strip() for l in line.strip().split('='))
                            if p[0] in ('name', 'identity', 'id'):
                                params['identity'] = p[1]
                            elif p[0] in ('color', 'colour'):
                                params['color'] = p[1]
                            elif p[0] in ('mass', 'radius'):
                                params[p[0]] = mpf(p[1])
                            elif p[0] == 'bounce':
                                params[p[0]] = float(p[1])
                            elif p[0] == 'position':
                                tp = p[1].strip('()[]\'').split(',')
                                params['init_pos'] = list(mpf(t) for t in tp)
                            elif p[0] == 'velocity':
                                tp = p[1].strip('()[]\'').split(',')
                                params['init_vel'] = list(mpf(t) for t in tp)
                            elif '[' in p[0] or ']' in p[0] or '*' in p[0] or '#' in p[0]:
                                pass
                            else:
                                tqdm.write(f'«obj_from()» → [!] Couldn\'t parse "{line}" on line {i} so it has been skipped.')
                        print(params)
                        bodies[key] = Body(**params)
                    # init engine and add bodies
                    engine = Engine(**eng_params)
                    engine.attach_bodies(list(bodies.values()))
                    # parse attributes
                    for var in var_strings:
                        if var[0] in ('do_collisions','do_bodygravity','do_fieldgravity'):
                            setattr(engine, val[0], eval(val[1]))
                        else:
                            tqdm.write(f'«obj_from()» → [!] Attempted to alter invalid attribute: {val[0]}.')
                    # parse callables
                    for val in func_strings:
                        if val[0] == 'plane' and len(val) == 3:
                            engine.create_plane(str(val[1]), float(val[2]))
                        elif val[0] == 'field' and len(parts) == 2:
                            vect = [mpf(v) for v in val[1].strip("'").replace('(', '').replace(')', '').split(',')]
                            engine.create_acceleration(vect)
                        elif val[0] == 'rel' and len(val) == 2:
                            engine.make_relative_to(bodies[val[1]])
                        elif val[0] == 'sim' and len(val) == 2:
                                engine.simulate(int(val[1]))
                        else:
                            tqdm.write(f'«obj_from()» → [!] Couldn\'t parse function "{' '.join(val)}" so it has been skipped.')
                    return engine
                elif obj == 'body':
                    # get body info and init
                    params = {}
                    for line in _input:
                        p = list(l.strip() for l in line.strip().split('='))
                        if p[0] in ('name', 'identity', 'id'):
                            params['identity'] = p[1]
                        elif p[0] in ('color', 'colour'):
                            params['color'] = p[1]
                        elif p[0] in ('mass', 'radius'):
                            params[p[0]] = mpf(p[1])
                        elif p[0] == 'bounce':
                            params[p[0]] = float(p[1])
                        elif p[0] == 'position':
                            tp = p[1].strip('()[]\'').split(',')
                            params['init_pos'] = list(mpf(t) for t in tp)
                        elif p[0] == 'velocity':
                            tp = p[1].strip('()[]\'').split(',')
                            params['init_vel'] = list(mpf(t) for t in tp)
                        else:
                            tqdm.write(f'«obj_from()» → [!] Couldn\'t parse "{line}" so it has been skipped.')
                    return Body(**params)
                elif obj == 'bodyext':
                    parsing_position = False
                    found_dt = False
                    positions = []
                    dt = 0
                    result = {}
                    # get body info
                    for line in _input:
                        parts = line.split('=')
                        if len(parts) == 2:
                            key, value = map(str.strip, parts)
                            value = value.replace('(', '').replace(')', '').replace(',', '')
                            if key in ('name', 'identity', 'id'):
                                result['identity'] = value
                            elif key in ('color', 'colour'):
                                result['color'] = (value if value is not 'None' else None)
                            elif key in ('mass', 'radius'):
                                result[key] = mpf(value)
                            elif key == 'bounce':
                                result[key] = float(value)
                            elif key == 'position':
                                parsing_position = True
                                continue
                            elif key == 'dt':
                                found_dt = True
                                dt = float(value)
                            elif key == 'velocity':
                                pass
                            else:
                                tqdm.write(f'«obj_from()» → [!] Couldn\'t parse "{line}" so it has been skipped.')
                        # check if line is start of position, gets pos data
                        if parsing_position == True:
                            if line.strip() == ')':
                                result['position'] = positions
                                parsing_position = False
                            else:
                                position_str = line.strip().replace('(', '').replace(')', '').replace("'", "")
                                position_tuple = tuple(map(mpf, position_str.split(',')))
                                positions.append(position_tuple)
                    # if dt specified, interpolate vel.
                    if found_dt == True:
                        params = {}
                        params['init_pos'] = result['position'][0]
                        posrec = result['position']
                        result.pop('position')
                        params.update(**result)
                        body = Body(**params)
                        for i,x in enumerate(posrec[1:]):
                            body.pos.next(x)
                            body.vel.next((Vector(posrec[i])-x)/dt)
                            body.acc.next((0,0,0))
                    return body
        else:
            e.raise_fnf_error(object)
    else:
        e.raise_type_error('object', str, object)

def export_obj(object, loc, overwritefile=True, info_name=None, csvs=True, direxists=True):
    # make file directory
    if isinstance(object,(Engine, Body)):
        os.makedirs(loc, exist_ok=direxists)
    if isinstance(object, Engine):
        if info_name == None:
            fname = 'eng_info'
        else:
            # create eng specific info string
            fname = info_name
        eng_info = [
            f'dt = {object.dt}\n',
            f'checking_range = {object._rangechk}\n',
            f'# number of intervals: {len(object)}\n',
            f'# bodies:({len(object.bodies)})\n',
            f'~ do_collisions = {object.do_collisions}\n',
            f'~ do_bodygravity = {object.do_bodygravity}\n',
            f'~ do_fieldgravity = {object.do_fieldgravity}\n'
        ]
        ids = []
        # append info for each object
        for b in object.bodies:  
            obj_id = str(b.identity).replace(' ','_').split(':')[0].split('(Static)')[0]
            ids.append(obj_id)
            for line in [f'*{obj_id} [\n',
                         f'name = {b.identity}\n',
                         f'mass = {b.mass.c()}\n',
                         f'radius = {b.radius.c()}\n',
                         f'color = {b.color}\n',
                         f'bounce = {b.bounce}\n',
                         f'position = {tup(b.pos)}\n',
                        f'velocity = {tup(b.vel)}\n' 
                         ']\n']:
                eng_info.append(line)
        # add planes and fields
        for p in object.planes:
            eng_info.append(f'-* plane = {p[0]} = {p[1]}\n')
        for v in object.fields:
            eng_info.append(f'-* field = {tup(v)}\n')
        # write to file
        with open(f'{loc}/{fname}.txt', ('w' if overwritefile is True else 'x')) as eng:
            eng.writelines(eng_info)
        # write csvs for each body
        if csvs == True:
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
                    d_df.to_csv(f'{loc}/{ids[i]}.csv', sep=',',index=False, mode=('w' if overwritefile is True else 'x'))
    if isinstance(object, Body):
        if info_name == None:
            fname = 'body_info'
        else:
            fname = info_name
        # make body info string
        obj_id = str(object.identity).replace(' ','_').split(':')[0].split('(Static)')[0]
        bod_info = [f'*{obj_id} [\n',
                    f'name = {object.identity}\n',
                    f'mass = {object.mass.c()}\n',
                    f'radius = {object.radius.c()}\n',
                    f'color = {object.color}\n',
                    f'bounce = {object.bounce}\n',
                    f'position = {tup(object.pos)}\n',
                    f'velocity = {tup(object.vel)}\n' 
                    ']\n']
        # write to file
        with open(f'{loc}/{fname}.txt', ('w' if overwritefile is True else 'x')) as bod:
            bod.writelines(bod_info)
        # make csv for body
        if csvs == True:
            info_dict ={'posX':object.pos.X.record,
                        'posY':object.pos.Y.record,
                        'posZ':object.pos.Z.record,
                        'velX':object.vel.X.record,
                        'velY':object.vel.Y.record,
                        'velZ':object.vel.Z.record,
                        'accX':object.acc.X.record,
                        'accY':object.acc.Y.record,
                        'accZ':object.acc.Z.record
                        }
            d_df = pd.DataFrame.from_dict(info_dict, 'index').transpose()
            d_df.to_csv(f'{loc}/{fname}.csv', sep=',',index=False, mode=('w' if overwritefile is True else 'x'))


