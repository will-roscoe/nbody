# NBody __init__.py
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
     |   ├── errors.py 
     |   ├── formatter.py 
     |   ├── horizons.py 
     |   └── io.py
     |
     ├── __init__.py
     └── config.py
   

Valid Import Sequence
     .errors
     .base
     .core.body
     .core.engine
     .tools.formatter
     .core.visual
     .tools.io
     .tools.horizons
"""
from .core.body import Body
from .core.engine import Engine
from .core.visual import MPLVisual
from .tools.horizons import horizons_batch, horizons_query
from .tools.io import obj_from, export_obj
