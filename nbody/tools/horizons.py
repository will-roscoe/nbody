import re
from astroquery.jplhorizons import Horizons
import astropy.units as ast
from datetime import datetime, timedelta

from tqdm import tqdm
from . import errors as e

def horizons_query(searchquery,observer='0',time='2023-11-03',num_type=float,return_type='body'):
    from ..core.body import typecheck, Body
    typecheck(((searchquery,str),(observer,str),(time,str),(num_type,type),(return_type,str)))
    if all([return_type.lower() != opt for opt in ('body','dict','print')]) :
        e.raise_value_error('return_type',type,return_type)
    
    tqdm.write(f'Querying "{searchquery}" @JPL Horizons System')
    _later_time = (datetime.strptime(time,'%Y-%m-%d')+timedelta(days=2)).strftime('%Y-%m-%d')
    
    _raw0 = Horizons(id=searchquery,location=observer,epochs={'start':time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors(get_raw_response=True)
    _tab = Horizons(id=searchquery,location=observer,epochs={'start':time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors()
    
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _raw0)
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        tqdm.write(f'could not find name for object "{searchquery}", reverting to placeholder name.')
        name = f'JPL_{searchquery}'
    
    try:
        _raw =_raw0.split('Ephemeris')[0]
        m_exp_unit = _raw.split('Mass')[1].split('^')[1].split('=')[0].strip(') ,~ (')
        m_exp = "".join(c for c in m_exp_unit if c.isdigit() or c == '.')
        m_unit = "".join(c for c in m_exp_unit.lower() if c == 'k' or c == 'g')
        mass = "".join(c for c in  _raw.split('Mass')[1].split('^')[1].split('=')[1].strip(') ,~ (').lower() if 
                        c.isdigit() or c == '.' or c=='+' or c=='-')
    except IndexError:
        print(_raw0)
    
    try:
        rad_string = _raw.lower().split('vol. mean radius')[1].split('=')[0:2]
    except IndexError:
        rad_string = _raw.lower().split('radius')[1].split('=')[0:2] 
    
    r_unit = "".join(c for c in rad_string[0].lower() if c == 'k' or c == 'm') 
    rad = list(c for c in rad_string[1].split(' ') if c != '') 
    _rad = []
    
    for x in rad:
        if any(char.isdigit() for char in x): 
            _rad.append(x)
    mass = num_type(''.join(mass.split('+-')[0]))*num_type(10**int(m_exp))
    rad = num_type(_rad[0].split('+-')[0])
    
    if r_unit == 'km':
        rad *= 1000
    if m_unit == 'g':
        mass /= 1000
    
    x, y, z = [num_type(_tab[pos].quantity[0].to_value(ast.m)) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [num_type(_tab[pos].quantity[0].to_value(ast.m/ast.s)) for pos in ('vx', 'vy', 'vz')]
    
    if return_type.lower() == 'print':
        print(f'''Object: {name}\n***********\nMass: {mass} kg\nRadius: {rad} m\n
***********\nposition: ({x}, {y}, {z}) (m)\nvelocity: ({vx}, {vy}, {vz}) (ms^-1)\n
***********\nQuery Date: {time}''')
    
    elif return_type.lower() == 'dict':
        return {'identity': name, 'mass': mass, 'radius': rad, 'position':(x,y,z), 'velocity':(vx,vy,vz)}
    
    elif return_type.lower() == 'body':
        return Body(mass=mass, init_pos=(x,y,z), init_vel=(vx,vy,vz), radius=rad, identity=name)


def horizons_batch(search_queries,observer='0',time='2023-11-03',num_type=float,return_type='body'):
    new_bodies = []
    for query in tqdm(search_queries,desc='Getting data from JPL Horizons',unit='queries'):
        if isinstance(query,str):
            new_bodies.append(horizons_query(query,observer,time,num_type,return_type))
        else:
            raise e.raise_type_error('item in search_queries',str,query)
    return new_bodies