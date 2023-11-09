from datetime import datetime, timedelta
import re

from tqdm import tqdm #pip install tqdm
from astroquery.jplhorizons import Horizons #pip install astroquery
import astropy.units as u # dependency of ^^^^ so should get installed then. 

#---------------------
# I have commented out functions/ classes related to other bits of code and replaced them with placeholders, 
# return_type = 'body' will not return anything due to this and errors are not descrptive.
#---------------------




def horizons_query(searchquery: str,
                   observer: str = '0',
                   time: str = '2023-11-03',
                   num_type: type = float,
                   return_type: str = 'dict'):
    """
    Create a Horizons object representing a celestial body using Horizons query data.

    This function queries NASA's Horizons system for information about a celestial body and creates a
    `Body` object with the retrieved data, including mass, initial position, initial velocity, radius,
    and identity.

    Args:
        searchquery (str): The identifier or name of the celestial body to query.
        observer (str, optional): The observer location (default is '0' for the solar system barycenter).
        time (str, optional): The date for which to retrieve data in the 'YYYY-MM-DD' format
                                (default is '2023-11-03').
        num_type (type, optional): Numeric Type to format values from source as. (default is float)
        return_type (str, optional): what format to return the object in. either 'print' to print output, 
                                    'dict' to return dictionary containing values or 'body' to return a Body
                                    instance. (default is 'body')
    Returns:
        Body: A `Body` object representing the celestial body with mass, initial position, initial velocity,
        radius, and identity information.
    OR  dict: dictionary with keys ('identity', 'mass', 'radius', 'position', 'velocity'). last 2 are 3-tuples.
    OR  print: prints a string representation of the dict above in a user readable form.
    Raises:
        TypeError: If one or more args are incorrect type.
        ValueError: if 'return_type' is not a valid option. make sure it is any of 'print', 'dict', 'body'.

    """
    # checking input argument types and/or values
    _args = ((searchquery, str),(observer,str), (time, str), (num_type, type), (return_type,str))
    for (inpt, typ) in _args:
        if not isinstance(inpt, typ):
            raise TypeError #e.raise_type_error('input', typ, inpt)
    if all([return_type.lower() != opt for opt in ('body', 'dict', 'print')]) :
        raise ValueError#e.raise_value_error('return_type', type, return_type)
    
    tqdm.write(f'Querying "{searchquery}" @JPL Horizons System')
        #find an arbitrary date after time that has chosen for 'epochs' dict
    _later_time = (datetime.strptime(time, '%Y-%m-%d') + timedelta(days=2)).strftime('%Y-%m-%d')
        #get raw data (for name, radius, mass) from JPL HORIZONS
    _raw = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors(get_raw_response=True)
        #get tabulated data in the form of an astropy Table from JPL HORIZONS
    _tab = Horizons(id=searchquery, location=observer, epochs={'start': time,
                                                               'stop':_later_time,
                                                               'step':'3d'}).vectors()
        # get name and send to a final return variable
    _n = re.search(r"Target body name: (.+?) \((\d+)\)", _raw)
    if _n:    
        name = f'{_n.groups()[0]} ({_n.groups()[1]})'
    else:
        tqdm.write(f'could not find name for object "{searchquery}", reverting to placeholder name.')
        name = f'JPL_{searchquery}'
    
    _raw =_raw.split('Ephemeris')[0] #remove useless data 

    # splits raw to get a string which will contain the mass exponent and unit, 
    # eg 10^24, kg and remove spaces and other punctuation
    m_exp_unit = _raw.split('Mass')[1].split('^')[1].split('=')[0].strip(') ,~ (')
    #pull out exponent from above
    m_exp = "".join(c for c in m_exp_unit if c.isdigit() or c == '.')
    # pull out unit from above
    m_unit = "".join(c for c in m_exp_unit.lower() if c == 'k' or c == 'g')
    #pull out mass from raw by splitting on =, will contain the value and whatever 
    # the next variable alongs' name is, remove other punctuation
    mass = "".join(c for c in  _raw.split('Mass')[1].split('^')[1].split('=')[1].strip(') ,~ (').lower() if 
                   c.isdigit() or c == '.' or c=='+' or c=='-')
    
    try: # try all known phrases that can represent radius, and get the unit and value
        rad_string = _raw.lower().split('vol. mean radius')[1].split('=')[0:2]
    except IndexError:
        # may pull multiple values out.
        rad_string = _raw.lower().split('radius')[1].split('=')[0:2] 
    
    # pull unit out
    r_unit = "".join(c for c in rad_string[0].lower() if c == 'k' or c == 'm') 
    #split on spaces, remove spaces from each item and remove empty items.
    rad = list(c for c in rad_string[1].split(' ') if c != '') 

    _rad = []
    for x in rad:
        #pull out only items that contain a number, will give a singular item which would
        #look something like '1234.0123+-0.067'
        if any(char.isdigit() for char in x): 
            _rad.append(x)
    
    # compute final return variables from parts
    mass = num_type(''.join(mass.split('+-')[0]))*num_type(10**int(m_exp))
    rad = num_type(_rad[0].split('+-')[0])
    # check if unit wasnt standard form, change finals if needed
    if r_unit == 'km':
        rad *= 1000
    if m_unit == 'g':
        mass /= 1000
    # get x,y,z position and velocity from tab data and convert to standard form
    x, y, z = [num_type(_tab[pos].quantity[0].to_value(u.m)) for pos in ('x', 'y', 'z')]
    vx, vy, vz = [num_type(_tab[pos].quantity[0].to_value(u.m/u.s)) for pos in ('vx', 'vy', 'vz')]
    
    # return options: dict, print a string, Body (custom class)
    
    if return_type.lower() == 'print': # print output
        print(f'''Object: {name}\n***********\nMass: {mass} kg\nRadius: {rad} m\n
***********\nposition: ({x}, {y}, {z}) (m)\nvelocity: ({vx}, {vy}, {vz}) (ms^-1)\n
***********\nQuery Date: {time}''')
    elif return_type.lower() == 'dict': # return a dictionary
        return {'identity': name, 'mass': mass, 'radius': rad, 'position':(x,y,z), 'velocity':(vx,vy,vz)}
    elif return_type.lower() == 'body': # return a Body instance
        pass #return Body(mass=mass, init_pos=(x,y,z), init_vel=(vx,vy,vz), radius=rad, identity=name)
    # no need to catch other possibilities as we have already checked at start