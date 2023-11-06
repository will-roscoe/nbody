from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta
import re
import astropy.units as u


def horizons_object(searchquery, observer='0', time='2023-11-03'):
    """
    Create a Horizons object representing a celestial body using Horizons query data.

    This function queries NASA's Horizons system for information about a celestial body and creates a
    `Body` object with the retrieved data, including mass, initial position, initial velocity, radius,
    and identity.

    Args:
        searchquery (str): The identifier or name of the celestial body to query.
        observer (str, optional): The observer location (default is '0' for the solar system barycenter).
        time (str, optional): The date for which to retrieve data in the 'YYYY-MM-DD' format (default is '2023-11-03').

    Returns:
        Body: A `Body` object representing the celestial body with mass, initial position, initial velocity,
        radius, and identity information.

    Raises:
        LookupError: If mass or radius information could not be found in the query output.
    """
    
    'Vol. Mean Radius (km) =  2440+-1'
    ''
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
    _raw = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors(get_raw_response=True)
    _tab = Horizons(id=searchquery, location=observer, epochs={'start': time, 'stop':_later_time,'step':'3d'}).vectors()
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _raw)
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        print('could not find name for object, reverting to placeholder name.')
        name = None
    _raw =_raw.split('Ephemeris')[0]
    m_exp_unit = _raw.split('Mass')[1].split('^')[1].split('=')[0].strip(') ,~ (')
    m_exp = "".join(c for c in m_exp_unit if c.isdigit() or c == '.')
    m_unit = "".join(c for c in m_exp_unit.lower() if c == 'k' or c == 'g')
    mass = "".join(c for c in  _raw.split('Mass')[1].split('^')[1].split('=')[1].strip(') ,~ (').lower() if c.isdigit() or c == '.' or c=='+' or c=='-')
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
    print(rad, mass, m_exp)
    rad = float(_rad[0].split('+-')[0])
    mass = float(mass[0].split('+-')[0])*10.**(float(m_exp))
    print(rad, mass)
    if r_unit == 'km':
        rad *= 1000
    if m_unit == 'g':
        mass /= 1000
    print(rad, mass)
    x, y, z = [_tab[pos].quantity[0].to_value(u.m) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [_tab[pos].quantity[0].to_value(u.m/u.s) for pos in ('vx', 'vy', 'vz')]
   
   
[horizons_object(x) for x in ('10', '199', '299', '499', '599', '699', '799', '899')]

