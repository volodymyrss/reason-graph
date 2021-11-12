#!/bin/env python

import os
import rdflib
import requests
import odakb.sparql
import time as time_module

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.io import fits

# def isgri_psla_arcmin(snr):
#     #Scaringi 2010
#     a=31.
#     b=0.25
#     c=-1.25

#     return a*snr**c+b


# class Catalogs(object):
#     simbad_disabled=False

#     def __init__(self,imcat=None,refcat_fn=None):
#         Simbad.remove_votable_fields("otype")
#         Simbad.add_votable_fields("otype")
#         Simbad.remove_votable_fields("coo_wavelength")
#         Simbad.add_votable_fields("coo_wavelength")

#         self.ISM_otype = "PoC PN? CGb bub EmO Cld GNe BNe DNe RNe MoC glb cor SFR HVC HII PN sh SR? SNR cir of? out HH"
#         self.ISM_otype += "ISM PartofCloud PN? ComGlob Bubble EmObj Cloud GalNeb BrNeb DkNeb RfNeb MolCld Globule denseCore SFregion HVCld HII PN HIshell SNR? SNR Circumstellar outflow? Outflow HH"
#         self.Star_otype = "* *iC *iN *iA *i* V*? Pe* HB* Y*O Ae* Em* Be* BS* RG* AB* C* S* sg* s*r s*y s*b HS* pA* WD* ZZ* LM* BD* N* OH* CH* pr* TT* WR* PM* HV* V* Ir* Or* RI* Er* Fl* FU* RC* RC? Ro* a2* Psr BY* RS* Pu* RR* Ce* dS* RV* WV* bC* cC* gD* SX* LP* Mi* sr* SN* su* Pl? Pl"
#         self.Star_otype += "Star *inCl *inNeb *inAssoc *in** V*? Pec* HB* YSO Ae* Em* Be* BlueStraggler RGB* AGB* C* S* SG* RedSG* YellowSG* BlueSG* HotSubdwarf post-AGB* WD* pulsWD* low-mass* brownD* Neutron* OH/IR CH pMS* TTau* WR* PM* HV* V* Irregular_V* Orion_V* Rapid_Irreg_V* Eruptive* Flare* FUOr Erupt*RCrB RCrB_Candidate RotV* RotV*alf2CVn Pulsar BYDra RSCVn PulsV* RRLyr Cepheid PulsV*delSct PulsV*RVTau PulsV*WVir PulsV*bCep deltaCep gammaDor pulsV*SX LPV* Mira semi-regV* SN Sub-stellar Planet? Planet"

#         self.imcat=imcat

#         self.load_refcat(refcat_fn)

#     def load_refcat(self,refcat_fn):
#         if refcat_fn is None:
#             refcat_fn = os.environ['ISDC_REF_CAT'][:-3]

#         self.refcat=fits.open(refcat_fn)[1].data
#         self.refcat_c=fits.open(refcat_fn)[1].data
#         self.refcat_coords = SkyCoord(self.refcat['RA_OBJ'], self.refcat['DEC_OBJ'], unit=u.deg)
#         self.refcat_coords = SkyCoord(self.refcat['RA_OBJ'], self.refcat['DEC_OBJ'], unit=u.deg)
#         self.refcat_names = self.refcat['NAME']

#         # self.refcat_df=helpers.fits_table_to_pandas(refcat_fn)

#         # self.refcat_df.NAME = np.array([kk.strip() for kk in self.refcat_df.NAME])#  map(str.strip, self.refcat_df.NAME))
    
#     def query_coord(self,coord,snr=5):
#         for i in reversed(range(10)):
#             try:
#                 return self._query_coord(coord, snr)
#             except requests.exceptions.ConnectionError as e:
#                 if i>0:
#                     print("retrying simbad!", i)
#                     time_module.sleep(1+i)
#                 else:
#                     raise

#     def query_name(self, name):

#     def _query_coord(self,coord,snr=5):
#         # for source,snr in sorted(zip(sources,source_snr),key=lambda x:-x[1]):
#         source = coord

#         rr = source.match_to_catalog_sky(self.refcat_coords)
#         if rr[1].deg < self.refcat_df.ERR_RAD[int(rr[0])] + isgri_psla_arcmin(snr) / 60.:
#             self.refcat_name = self.refcat_names[rr[0]]
#         else:
#             self.refcat_name = "unmatched"

#         self.all_associations = []
#         self.all_otypes = []
#         self.interesting_associations = []
#         self.interesting_otypes = []
#         self.integral_cat_source = False
#         self.best_names = []
#         self.best_names_all = []
#         self.aka = []

#         if self.simbad_disabled:
#             return

#         r=None
#         for ntry in range(5):
#             try:
#                 r = Simbad.query_region(source, radius=getattr(self,'search_radius_arcmin',10.) * u.arcmin)
#                 break
#             except Exception as e:
#                 print("simbad refused!",source,e)
#                 pass
#             time_module.sleep(3)

#         if r is None:
#             print("simbad returns nothing")
#             self.imcat_name = "unmatched"
#             return
#         else:
#             print("simbad returns",len(r))

#         #???rxi = Simbad.query_objectids(r['MAIN_ID'])

#         #  print "extra ids:",len(rxi)

#         #???rx = Simbad.query_object(rxi['ID'])
#         # r=astropy.table.vstack([r,])
#         # print set(rx['OTYPE'])

#         # display(r)

#         r['distance'] = [round(x, 2) for x in
#                          source.separation(SkyCoord(r['RA'], r['DEC'], unit=(u.hourangle, u.deg))).arcmin]

#         r.remove_columns(
#             ["RA_PREC", "DEC_PREC", "COO_ERR_MAJA", "COO_ERR_MINA", "COO_QUAL", "COO_BIBCODE", "COO_WAVELENGTH",
#              "COO_WAVELENGTH_2"])

#         m_uninteresting = (r['OTYPE'] == "ISM") | \
#                           (r['OTYPE'] == "UV") | \
#                           (r['OTYPE'] == "Star") | \
#                           (r['OTYPE'] == "Galaxy") | \
#                           (r['OTYPE'] == "PM*") | \
#                           (r['OTYPE'] == "DkNeb") | \
#                           (r['OTYPE'] == "RRLyr") | \
#                           (r['OTYPE'] == "**") | \
#                           (r['OTYPE'] == "YSO") | \
#                           (r['OTYPE'] == "Candidate_YSO") | \
#                           (r['OTYPE'] == "Bubble") | \
#                           (r['OTYPE'] == "Radio") | \
#                           (r['OTYPE'] == "HI") | \
#                           (r['OTYPE'] == "Candidate_AGB*") | \
#                           (r['OTYPE'] == "Maser") | \
#                           (r['OTYPE'] == "Radio(mm)") | \
#                           (r['OTYPE'] == "Radio(sub-mm)") | \
#                           (r['OTYPE'] == "Mira") | \
#                           (r['OTYPE'] == "semi-regV*") | \
#                           (r['OTYPE'] == "Cl*") | \
#                           (r['OTYPE'] == "IR")

#         for ot in set(r[~m_uninteresting]['OTYPE']):
#             if any([ot == o for o in
#                     self.ISM_otype.split() + self.Star_otype.split()
#                     ]) or \
#                     ("*" in ot and ot != "Symbiotic*" and ot != "CV*"):
#                 m_uninteresting |= (r['OTYPE'] == ot)

#         self.all_associations = r
#         self.all_otypes = list(set(r['OTYPE']))
#         self.interesting_associations = r[~m_uninteresting]
#         self.interesting_otypes = list(set(r['OTYPE'][~m_uninteresting]))

#         self.integral_cat_source = False
#         self.best_names = []
#         self.aka = []

#         def pick_names(rr):
#             if len(rr) > 0:
#                 if len(rr['MAIN_ID']) < getattr(self,'max_query',100):
#                     aka = list(Simbad.query_objectids(rr['MAIN_ID'])['ID'])

#                     favs = ["Granat", "GX", "GRO", "IGR", "GRS", "4U", "XTE", "Ginga", "1A", "^3A", "OAO", "^1E"]
#                     disfavs = ["CC", "GRS G"]
#                     best_names = [ak for ak in aka if any([re.search(fav, ak) for fav in favs]) and not any(
#                         [disfav in ak for disfav in disfavs])]
#                     return best_names,aka
#                 else:
#                     print("refusing to query: too large",len(rr['MAIN_ID']))
#             return [],[]


#         self.best_names, self.aka = pick_names(r[~m_uninteresting])
#         self.best_names_all, self.aka_all = pick_names(r)
#         self.integral_cat_source = any(["INT" in ak for ak in self.aka])


#         if self.imcat is not None:
#             ri = source.match_to_catalog_sky(self.imcat_coords)
#             if ri[1].arcmin < 30:
#                 self.imcat_name = self.imcat_names[float(ri[0])]
#             else:
#                 self.imcat_name = "unmatched"
#         else:
#             self.imcat_name=""

#     # def display(self):
#     #     display(HTML("""<h1><font style=\"color: #aa0000; font-size:20px\"> %s</font></h1>
#     #                     <h2><font style=\"color: #00aa00\"> %s</font>
#     #                         <font style=\"color: #ee0000\"> %s</font></h2>
#     #                     <h2><font style=\"color: #0000aa\"> %s</font></h2>
#     #                         """%(
#     #                      "Image: " + self.imcat_name + ", REF: " + self.refcat_name + ", " +
#     #                      (self.best_names[0] if len(self.best_names) > 0 else "?"),
#     #                      "INTEGRAL known" if self.integral_cat_source else "",
#     #                      "INTEGRAL UNKNOWN!" if not self.integral_cat_source and self.refcat_name == "unmatched" else "",
#     #                      "/".join(self.interesting_otypes),
#     #                  )))

#     #     display(self.interesting_associations)

#     #     display(HTML("<h3 style=\"color: #337766\"> AKA %s </h3>"
#     #                  % (
#     #                      ", ".join(self.best_names),
#     #                  )))

#     #     display(pd.DataFrame([self.aka[i * 5:i * 5 + 5] for i in range(len(self.aka) / 5 + 1)]))


# Catalogs()

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

    # print(r)

    if 'EXO' not in w['label']:
        continue
    
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


