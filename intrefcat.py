import re
from astropy.io import fits as fits
import os

import numpy as np 
os.system('wget -c ftp://isdcarc.unige.ch/arc/rev_3/cat/hec/gnrl_refr_cat_0043.fits')

cat = fits.open("gnrl_refr_cat_0043.fits")[1].data

m = cat['ISGRI_FLAG2'] == 5

import rdflib
import odakb.sparql

G = rdflib.Graph()
isdcrefcat_ns = rdflib.Namespace('http://odahub.io/ontology/intrefcat#')
oda_ns = rdflib.Namespace('http://odahub.io/ontology#')
rdfs_ns = rdflib.Namespace('http://www.w3.org/2000/01/rdf-schema#')
rdf_ns = rdflib.Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')


# G.bind('oda', oda_ns)
# G.bind('isdcrefcat', isdcrefcat_ns)

for r in cat[m]:
    name = r['NAME']
    sname = re.sub("[^0-9a-zA-Z]", "", name)

    G.add((
        oda_ns[f'object/{sname}'],
        rdf_ns['type'],
        oda_ns['AstrophysicalObject']
    ))
    G.add((
        oda_ns[f'object/{sname}'],
        rdfs_ns['label'],
        rdflib.Literal(name)
    ))
    

    for c in cat.columns:
        if len(c.dtype.shape)>0:
            continue
        print(c.name, c.dtype, r[c.name])

        G.add((
            oda_ns[f'object/{sname}'],
            isdcrefcat_ns[c.name],
            rdflib.Literal(r[c.name])
        ))


S = G.serialize(format='turtle')
print(S)

odakb.sparql.insert("\n".join([
    f"{s.n3()} {p.n3()} {o.n3()} ." for s, p, o in G
]))