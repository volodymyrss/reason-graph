import odakb.sparql

for w in odakb.sparql.select('?wfl oda:isRequestingAstroObject ?obj'): 
    # this is not necessarily true, but will be applicable until we understand plan situation and make ontology of inputs
    r = f'''
            <{w['wfl']}> a oda:Workflow;
                         oda:has_input oda:AstroObject .
        '''

    print(r)

    odakb.sparql.insert(r)
