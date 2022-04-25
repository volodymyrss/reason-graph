import hashlib
import odakb.sparql
import re


def run(ingest=False):
    data = ""

    for event in odakb.sparql.select('?paper paper:mentions_named_event ?name', limit=10000): 
        if "PKS" not in event["name"]: continue

        print(event)

        r_uri = re.sub("[^0-9a-zA-Z]", "", event["name"])

        r = f'''
            oda:{r_uri} a oda:AstrophysicalObject;
                              rdfs:label "{event["name"]}" .

            <{event['paper']}> paper:mentions_named_event oda:{r_uri} .
       
            oda:{r_uri} oda:cited_in <{event['paper']}>  .
            '''

        print(r)

        data += r

    odakb.sparql.insert(data)


if __name__ == "__main__":
    run(ingest=True)
