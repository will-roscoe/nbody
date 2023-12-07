def raise_type_error(obj:str, typ: tuple|object, invalid_obj:object):
    '''raises the TypeError:
{invalid_obj} is not a valid {typ} instance for {str(obj)}, type given was {type(invalid_obj)}.'''
    raise TypeError(f'{invalid_obj} is not a valid {typ} instance for {str(obj)},\
type given was {type(invalid_obj)}.')

def raise_len_error(obj:str, correct_len:int, invalid_obj:object):
    '''raises the IndexError: {invalid_obj} of length {len(invalid_obj)} was given for {str(obj)},
which should have length of {correct_len}.'''
    raise IndexError(f'{invalid_obj} of length {len(invalid_obj)} was given for {str(obj)}, \
which should have length of {correct_len}.')

def raise_value_error(obj:str, typ: tuple|object, invalid_obj:object):
    '''raises the ValueError:
{obj} of type {typ} was passed a invalid value of {invalid_obj}.'''
    raise ValueError(f'{obj} of type {typ} was passed a invalid value of {invalid_obj}.')

def raise_component_error(obj:str, invalid_obj:tuple|list):
    '''raises the IndexError:
input tuple for {obj} does not have 3 components, length was {len(invalid_obj)}.'''
    raise IndexError(f'input tuple for {obj} does not have 3 components, \
length was {len(invalid_obj)}.')

def raise_list_type_error(obj:str, typ: tuple|object, invalid_component:object):
    '''raises the TypeError: one or more iems in {str(obj)} is not a valid {typ}
instance for {str(obj)}, part raised was {invalid_component}.'''
    raise TypeError(f'one or more iems in {str(obj)} is not a valid {typ} \
instance for {str(obj)}, part raised was {invalid_component}.')

def raise_out_of_range(obj:str, arg):
    '''raises the IndexError: Index of {arg} is out of range for object {obj}.'''
    raise IndexError(f'Index of {arg} is out of range for object {obj}.')

def raise_read_only(obj):
    '''raises the IndexError:
{obj} is read and append only.'''
    raise IndexError(f'{obj} is read and append only.')

def raise_unmatched_error(objs, objs_insts):
    '''raises the RuntimeError: vector component lengths are unmatched:
[f'len{ostr}={len(oinst)}'for ostr, oinst in (objs, objs_insts)].'''
    temp = [f'len{ostr}={len(oinst)}'for ostr, oinst in (objs, objs_insts)]
    raise RuntimeError(f'vector component lengths are unmatched: {temp}.')

def raise_evaluation_error(objs):
    '''raises the RuntimeError: 'could not accurately evaluate body even after \
interpolating to give missing values.\
[f'{i}:{len(obj)}, {obj}.\n' for i, obj in enumerate(objs)].'''
    temp = [f'{i}:{len(obj)}, {obj}.\n' for i, obj in enumerate(objs)]
    raise RuntimeError(f'could not accurately evaluate body even after \
interpolating to give missing values. {temp}.')

def raise_fnf_error(loc):
    '''raises the FileNotFoundError: Could not find file at {loc}.'''
    raise FileNotFoundError(f'Could not find file at {loc}.')



