import hashlib
import odakb.sparql

for grb in odakb.sparql.select('?paper paper:mentions_named_grb ?name; paper:grb_isot ?grb_isot'): 
    time_instant = "oda:TimeI" + hashlib.md5(grb["grb_isot"].encode()).hexdigest()[:8]
    
    r = f'''
        oda:{grb["name"]} a oda:AstrophysicalObject;
                          rdfs:label "{grb["name"]}";
                          a oda:TransientObject;
                          oda:event_time_isot "{grb["grb_isot"]}";
                          oda:event_time {time_instant} .
        {time_instant} oda:isot "{grb["grb_isot"]}";
                       a oda:TimeInstant .
        '''

    print(r)

    odakb.sparql.insert(r)
