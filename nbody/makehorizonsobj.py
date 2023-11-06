
from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta
import astropy.units as u
import re
from body import Body
def make_horizons_object(searchquery, observer='0', time='2023-11-03'):
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
    _object = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors(get_raw_response=True)
    _object_tab = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors()
    _m = re.search(r"\s+?Mass\s+?([\S]+).*?=\s+?([\S]+?)\s+?", _object)
    _r = re.search(r'\s+?radius\s+?([\S]+).*?=\s+?([\S]+?)\s+?', _object)
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _object)
    if _m:  
        _mbase, _mexp_s = _m.groups()[1], _m.groups()[0]
        mass = float(_mbase.split('+-')[0]) * 10**(float(_mexp_s.split('x10^')[1]))
    else: 
        raise LookupError('could not find Mass in output')
    if _r:
        _rad =_r.groups()[1]
        radius = float(_rad.split('+-')[0])*1000
    else:
        raise LookupError('could not find radius (Vol. mean Radius) in output')
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        print('could not find name for object, reverting to placeholder name.')
        name = None
    x, y, z = [_object_tab[pos].quantity[0].to_value(u.m) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [_object_tab[pos].quantity[0].to_value(u.m/u.s) for pos in ('vx', 'vy', 'vz')]
   
    return Body(mass=mass, init_pos=(x,y,z), init_vel=(vx,vy,vz), radius=radius, identity=name)