# NBody __init__.py
'''
The goal of this project is to efficiently simulate the gravitational and kinetic interactions between physical objects
and display them in a user friendly interface. The usage of this project are object oriented and meant to be flexible 
and easy to use.

'''
"""
File Structure:
⸺ nbody/
     ├── core/
     |   ├── __init__.py 
     |   ├── base.py 
     |   ├── body.py 
     |   ├── engine.py 
     |   └── visual.py
     |   
     ├── examples/
     |   └── data/
     |  
     ├── tests/
     |   ├── testdata/
     |   |   ├── test_body.txt
     |   |   ├── test_bodyext.txt
     |   |   ├── test_engine.txt
     |   |   └── test_enginerun.txt
     |   ├── __init__.py 
     |   ├── test_base.py 
     |   ├── test_core.py 
     |   └── test_tools.py 
     |
     ├── tools/
     |   ├── __init__.py 
     |   ├── _config.py 
     |   ├── errors.py 
     |   ├── formatter.py 
     |   ├── horizons.py 
     |   └── io.py
     |
     ├── __init__.py
     ├── LICENSE.txt
     ├── MATHTYPE
     └── requirements.txt
"""
from .core.base import math_conf as CONFIG #noqa
from .core.body import Body
from .core.engine import Engine
from .core.visual import MPLVisual
from .tools.horizons import horizons_batch, horizons_query
from .tools.io import obj_from, export_obj
