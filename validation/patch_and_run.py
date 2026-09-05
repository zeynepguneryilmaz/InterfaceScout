from pathlib import Path
import runpy
src = Path('validation/run_method_validation.py').read_text()
src = src.replace("def topkeys(a,k=10):return [x[0] for x in sorted(a.items(),key=lambda z:(-z[1],z[0]))[:k] if xval(z:=z)>0]\ndef xval(z): return z[1]", "def topkeys(a,k=10):\n    return [item[0] for item in sorted(a.items(), key=lambda z:(-z[1],z[0])) if item[1] > 0][:k]")
tmp = Path('validation/_run_method_validation_fixed.py')
tmp.write_text(src)
runpy.run_path(str(tmp), run_name='__main__')
