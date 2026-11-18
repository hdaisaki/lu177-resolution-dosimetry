"""Provenance record of the figures (figure file -> generating script, panel count, md5)."""
import json, hashlib, os, datetime
from resdos.paths import ROOT
def record(fig, script, panels):
    path=os.path.join(ROOT,'results','figs','manifest.json')
    m=json.load(open(path)) if os.path.exists(path) else {}
    full=fig if os.path.isabs(fig) else os.path.join(ROOT,'results',fig)
    m[os.path.basename(fig)]=dict(script=script,panels=panels,md5=hashlib.md5(open(full,'rb').read()).hexdigest(),time=datetime.datetime.now().isoformat(timespec='seconds'))
    json.dump(m,open(path,'w'),indent=1)
