#!/bin/env python

import os
import rdflib
import requests
import odakb.sparql
import time as time_module

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.io import fits


from astroquery.simbad import Simbad

Simbad.remove_votable_fields("otype")
Simbad.add_votable_fields("otype")
Simbad.remove_votable_fields("coo_wavelength")
Simbad.add_votable_fields("coo_wavelength")


# refer simbad types to to https://ivoa.net/documents/Notes/AstrObjectOntology/
# https://www.ivoa.net/rdf/object-type/2020-10-06/object-type.rdf

# e.g. http://simbad.u-strasbg.fr/simbad/otypes#
# http://simbad.u-strasbg.fr/simbad/otypes#HXB

ivoa_ot_G = rdflib.Graph()

# why does it not work by url?
ivoa_ot_G.parse('https://www.ivoa.net/rdf/object-type/2020-10-06/object-type.rdf', format='xml')
#ivoa_ot_G.parse(requests.get('https://www.ivoa.net/rdf/object-type/2020-10-06/object-type.rdf').text, format='xml')


for w in odakb.sparql.select('?aobj a oda:AstrophysicalObject; rdfs:label ?label'): 
    # r = f'''
    #         <{w['wfl']}> a oda:Workflow;
    #                      oda:has_input oda:AstroObject .
    #     '''


    if 'Mrk' not in w['label']:
        continue

    print(f"\033[31m{w['aobj']}\033[0m")
    print(f"\033[31m{w['label']}\033[0m")
    
    result_table = Simbad.query_object(w['label'], wildcard=False)
    print(result_table)

    if result_table is None:
        continue

    otype = str(result_table[0]['OTYPE']).strip()

    print(f"\"{otype}\"")

    # TODO: deeper

    G = rdflib.Graph()
    for t in ivoa_ot_G.query(f'CONSTRUCT WHERE {{ ?obj rdfs:label "{otype}"; ?p ?o }}'):
        G.add(t)

    print(G.serialize(format='turtle'))

    objs = list(ivoa_ot_G.query(f'SELECT ?obj WHERE {{ ?obj rdfs:label "{otype}" }}'))
    
    r = f'''
            <{w['aobj']}> oda:simbadOTYPE "{otype}" .
        '''

    for obj in objs:
        r += f'''
            <{w['aobj']}> a <{obj[0]}> .
        '''

    for s, p, o in G:
        r += f'''
            {s.n3()} {p.n3()} {o.n3()} .
        '''



    print(r)

    odakb.sparql.insert(r)
    

    #break


