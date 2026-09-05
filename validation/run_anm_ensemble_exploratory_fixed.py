import json
import runpy

_original_dumps = json.dumps

def _safe_dumps(obj, *args, **kwargs):
    kwargs.setdefault('default', str)
    return _original_dumps(obj, *args, **kwargs)

json.dumps = _safe_dumps
runpy.run_path('validation/run_anm_ensemble_exploratory.py', run_name='__main__')
