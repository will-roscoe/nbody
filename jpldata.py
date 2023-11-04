from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u

object = Horizons(id='Sun', location='@0', epochs={'start':'2023-11-03', 'stop':'2023-11-04','step':'10d'})
tab = object.vectors()
_pos = [tab['x'].quantity, tab['y'].quantity, tab['z'].quantity]
_vel = [tab['vx'].quantity, tab['vy'].quantity, tab['vz'].quantity]
_mass = 



