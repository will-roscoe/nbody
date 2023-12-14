import pytest

import shutil
import os

from .. import export_obj, obj_from, Engine, Body, mplVisual
from ..tools.horizons import horizons_batch, horizons_query

def cleanup(folder):
    if os.path.exists(folder):
            shutil.rmtree(folder)
data_ = 'nbody/tests/testdata'
tmp_ = 'nbody/tmp/tools'
floc = {
    'body':'test_body.txt',
    'bodyext':'test_bodyext.txt',
    'eng':'test_engine.txt',
    'engrun':'test_enginerun.txt',
}
tloc = {
    'body':['tmp_body','.txt'],
    'bodyext':['tmp_bodyext','.txt'],
    'eng':['tmp_engine','.txt'],
}
queries = ('10','199','299','399','499','599','699','799','899')

class TestHorizonsAsDict:
    def test_get_single_dict(self):
        body = horizons_query(queries[1], return_type='dict')
        assert body['identity'] == 'Mercury (199)' and body['radius'] == 2440000.
    def test_get_batch_dict(self):
        bodies = horizons_batch(queries, return_type='dict')
        assert len(bodies) == 9


class TestHorizonsAsBody:
    def test_get_single_body(self):
        body = horizons_query(queries[1])
        assert isinstance(body, Body) and body.identity == 'Mercury (199)'
    def test_get_batch_body(self):
        bodies = horizons_batch(queries)
        assert all(isinstance(body, Body) for body in bodies) and len(bodies) == 9


class TestIO:
    def test_simple_body(self):
        try:
            body = obj_from(f'{data_}/{floc["body"]}', obj='body')
            _orig = body['_data_']
            export_obj(body, tmp_, True, tloc['body'][0], False)
            del body
            body = obj_from(f'{tmp_}/{''.join(tloc["body"])}', obj='body')
            _parsed = body['_data_']
        finally:
            cleanup(tmp_)
        assert all(_orig[key] == _parsed[key] for key in list(_orig.keys()))
    
    def test_bodyext(self):
        try:
            body = obj_from(f'{data_}/{floc["bodyext"]}', obj='bodyext')
            _orig = body['_data_']
            export_obj(body, tmp_, True, tloc['bodyext'][0], False)
            del body
            body = obj_from(f'{tmp_}/{''.join(tloc["bodyext"])}', obj='body')
            _parsed = body['_data_']
        finally:
            cleanup(tmp_)
        assert all(_orig[key] == _parsed[key] for key in list(_orig.keys()))
        
    def test_simple_engine(self):
        try:
            eng = obj_from(f'{data_}/{floc["eng"]}', obj='engine')
            _orig = eng['_data_']
            _or_b = [_ob['_data_'] for _ob in _orig.pop('bodies')]
            export_obj(eng, tmp_, True, tloc['eng'][0], False)
            del eng
            eng = obj_from(f'{tmp_}/{''.join(tloc["eng"])}', obj='engine')
            _parsed = eng['_data_']
            _pr_b = [_ob['_data_'] for _ob in _parsed.pop('bodies')]
        finally:
            cleanup(tmp_)
        _body_dup = all(all(_ob[key] == _pb[key] for key in list(_ob.keys())) for (_ob, _pb) in zip(_or_b, _pr_b))
        assert all(_orig[key] == _parsed[key] for key in list(_orig.keys())) and _body_dup

    def test_run_engine(self):
        eng = obj_from(f'{data_}/{floc["engrun"]}', obj='engine')
        assert len(eng) == 100