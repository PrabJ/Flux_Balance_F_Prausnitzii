# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # Integrated Bioinformatics Project : Flux_Balance_F_Prausnitzii 
# #### Students: Michael Shawn Neilsson, Biancamaria FLorenzi, Prabhat Juyal, Tim Blokker
# #### Supervisors: Clémence Joseph, Karoline Faust

# %%
#import packages
import io
import cobra
import functions
from functions import *


# %%
#model = cobra.io.read_sbml_model('Data/f_prausnitzii.xml')
#info(model)
model=cobra.io.load_matlab_model('Data/iFpraus_v_1_0.mat')
info(model)


# %%
outputmodel= cobra.io.save_json_model(model, "Data/iFpraus.json")


# %%
#unaltered model
model1=model.copy()
print("\nModel Medium")
#model1.objective = {model1.reactions.get_by_id('EX_ac(e)'): 1}
model1.objective = {model1.reactions.get_by_id('Biomass_FP'): 1}
print(model1.summary())

# %% [markdown]
# ### What are the different reactions
# e: extracellular
# 
# c: cytosol
# 
# EX_ : exchange reaction
# ### Gifu Anaerobic Medium (mGAM) from:
# https://hyserve.com/files/05433_GAM-Broth_Modified_final.pdf
# 
# 
# - Peptone, Soya Peptone, Proteose Peptone: source of amino acids -> most likely all of them 
# - Yeast extract, liver extract, meat extract, digested serum -> amino acids but also a lot of other crap, will leave this till the end 
# - Dextrose -> same as glucose (https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17634) and so -> glc_D
# - Soluble starch -> starch1200 (potatoe starch)
# - L-Tryptophane -> trp__L
# - L-Cystein Hydrochloride -> cys_L and cl and h (hydrochloride<-> HCL becomes cl and h)
# - Sodium Thioglycolate -> thiog and na1
# - L-Arginine -> arg_L
# - Vitamine K1 -> phllqne https://pubchem.ncbi.nlm.nih.gov/compound/Vitamin-K1#section=Depositor-Supplied-Synonyms 
# - Hemin -> Iron and Chlorine  fe3 and cl  and perhaps ppp9 https://pubchem.ncbi.nlm.nih.gov/#query=3-%5B18-(2-carboxyethyl)-8%2C13-bis(ethenyl)-3%2C7%2C12%2C17-tetramethylporphyrin-21%2C23-diid-2-yl%5Dpropanoic%20acid%3Biron(3%2B)%3Bchloride 
# - Potassium Dihydrogen Phosphate -> pi, k, h
# - Sodium Chloride -> na1 cl
# 
# Digested serum : https://www.nebiolabs.com.au/-/media/catalog/datacards-or-manuals/p8108datasheet-lot0021306.pdf	
# 
# Generally useful page ; https://pubchem.ncbi.nlm.nih.gov/#query=Hemin
# 
# check also: https://opencobra.github.io/cobratoolbox/latest/tutorials/tutorialMetabotoolsI.html 

# %%
import numpy as np
aa=("ala,arg,asn,asp,asx,cys,glu,gln,glx,gly,his,ile,leu,lys,met,phe,pro,ser,thr,trp,tyr,val").split(",")
aa=[aminoacid +"(e)" for aminoacid in aa]
print(aa)


# %%
MediumMgam = aa+['glc(e)', 'starch1200(e)', 'cl(e)', 'h(e)',"ac(e)", "so4(e)" ,'thiog(e)', 'na1(e)', 'fe3(e)', 'ppp9(e)', 'h2o(e)', 'pi(e)', 'k(e)']
print(MediumMgam)


# %%
model2=model.copy()
model_mgam=medium(model2,MediumMgam)
print('\nMedium_mgam')
for i in range(1,len(model_mgam.reactions)):
    model_mgam.reactions[i].upper_bound=1000
    if model_mgam.reactions[i].lower_bound == 0:
        model_mgam.reactions[i].lower_bound=0 #does not change till 77,78
    else:
        model_mgam.reactions[i].lower_bound=-1000 #does not change till 77,78
model_mgam.summary()


# %%
def find_match(medium_list):
    matching_id=[]
    not_matching_id=[]
    for component in medium_list:
        try:
            matching_id.append(model.exchanges.get_by_id("EX_"+component))
        except KeyError:
            not_matching_id.append(component)
            continue
    return matching_id, not_matching_id


# %%
find_match(MediumMgam)


# %%
fba = model.optimize()
fluxes = deleteNull(fba.fluxes)
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes),len(fluxes)/len(model.reactions)*100))
print(fba.objective_value)


# %%
from cobra.medium import minimal_medium

max_growth = model.slim_optimize()
minimal_medium(model, max_growth)


# %%
from cobra.medium import minimal_medium
max_growth = model_mgam.slim_optimize()
minimal_medium(model_mgam, max_growth)

# %% [markdown]
# ---
# ## Reinforced Clostridial Medium (RCM)
# http://www.oxoid.com/UK/blue/prod_detail/prod_detail.asp?pr=CM0149&org=53&c=UK&lang=EN
# 
# * Yeast Extract -> [lots of stuff in it, unsure of what's garbage](https://www.chemicalbook.com/ChemicalProductProperty_EN_CB9440339.htm#:~:text=Yeast%20extract%20has%20a%20protein,aromatic%20compounds%20and%20other%20components.)
#     * Glutathione: gthrd (reduced version)
#     * 18 AA's: see above
#     * Dextran: glc__D (it's made of glucose)
#     * Mannan: mannan
#     * Trehalose: tre
#     * B-vitamins: (Clemence mentioned vitamins being important so I'm putting them all in)
#         * \1. Thiamin -> thm
#         * \2. Riboflavin -> ribflv
#         * \3. Niacin -> trp__L (made of tryptophan)
#         * \5. Pantothenic Acid -> pnto__R
#         * \6. Pyridoxine -> pydxn
#         * \7. Biotin -> btn
#         * \9. Folic Acid -> ... can't find anything
#         * \12. Cobalamin -> b12 (or cbl1)
#     * Biotin
# * Peptone -> Sticking with the above logic, probably a source of all AAs
# * Glucose -> Assuming D_Glucose: glc__D
# * Soluble Starch -> starch1200
# * Sodium Chloride -> na1, cl
# * Sodium Acetate -> Assuming breaking down into sodium and acetate: ac, na1
# * Cysteine Hydrochloride -> cys_L, cl, h (could also be cys_D? cys_L is more widely used)
# * Agar: Agarose (70%) and Agaropectin (30%) [according to Wiki](https://en.wikipedia.org/wiki/Agar)
#     * Agarose: D-galactose & 3,6-anhydro-L-galactopyranose
#         * D-galactose: gal
#         * 3,6-anhydro-L-galactopyranose: [similar structure & compositon to beta-D-allose](https://pubchem.ncbi.nlm.nih.gov/#query=CID67020466%20structure&tab=similarity)
#             * D-Allose: all__D
#     * Agaropectin: D-galactose, L-galactose, pyruvate, sulfate
#         * D-galactose: gal
#         * L-galactose: gal__L
#         * Pyruvate: pyr
#         * Sulfate: so4
# %% [markdown]
# | Component              | Concentration (g/L) | BiGG Metabolites                                                                          |
# |------------------------|---------------------|-------------------------------------------------------------------------------------------|
# | Yeast Extract          | 13.0                | (see all AA's), gthrd, glc\__D, mannan, tre, thm, ribflv, trp\__L, pnto\__R, pydxn, btn, b12 |
# | Peptone                | 10.0                | (see all AA's)                                                                            |
# | Glucose                | 5.0                 | glc_D                                                                                     |
# | Soluble Starch         | 1.0                 | starch1200                                                                                |
# | Sodium Chloride        | 5.0                 | na1, cl                                                                                   |
# | Sodium Acetate         | 3.0                 | na1, ac                                                                                   |
# | Cysteine Hydrochloride | 0.5                 | cys_L, cl, h                                                                              |
# | Agar                   | 0.5                 | gal, all\__D, gal\__L, pyr, so4                                                           |

# %%
medium_rcm = list(('gthrd','glc__D','mannan','tre','thm','ribflv','trp__L','pnto__R','pydxn','btn','b12',
             'starch1200','na1','ac','cl','cys_L','gal','all__D','gal__L','pyr','so4'))
medium_rcm=[mets+"(e)" for mets in medium_rcm]
model_rcm=medium(model.copy(), medium_rcm)
for i in range(1,len(model_rcm.reactions)):
    model_rcm.reactions[i].upper_bound=1000
    model_rcm.reactions[i].lower_bound=-1000
model_rcm.objective = {model_rcm.reactions.get_by_id('Biomass_FP'): 1}
model_rcm.summary()


# %%
model_mmcb.reactions.get_by_id('mannan')


# %%
# Post-optimization
fba_rcm = model_rcm.optimize()
fluxes = deleteNull(fba_rcm.fluxes)
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes),len(fluxes)/len(model.reactions)*100))
print(model_rcm.objective.value)

# %% [markdown]
# ---
# ## mMCB
# 
# (There are also supplements mentioned in the paper but this is the base medium)
# * [Bacteriological Peptome (6.5) (Oxoid)](http://www.oxoid.com/UK/blue/prod_detail/prod_detail.asp?pr=LP0037&c=UK&lang=EN): 
#     * Vague, decided to go with polypeptides: <b>polypep</b>
#     * Nitrogen: <b>n2</b>
# * [Soy Peptome (5.0) (Oxoid)](http://www.oxoid.com/UK/blue/prod_detail/prod_detail.asp?pr=LP0044&cat=&c=UK&lang=EN) --> another vague one:
#     * Stachyose: <b>stys</b>
#     * Raffinose: <b>raffin</b>
#     * Sucrose: <b>sucr</b>
#     * Nitrogen: <b>n2</b>
# * [Yeast Extract (3.0) (VWR International, Darmstadt, Germany)](https://us.vwr.com/store/product/7437401/vwr-life-science-yeast-extract-bacteriological-grade):
#     * Vitamin B (same logic as RCM):
#         * 1. Thiamin -> <b>thm</b>
#         * 2. Riboflavin -> <b>ribflv</b>
#         * 3. Niacin -> <b>trp__L</b> (made of tryptophan)
#         * 5. Pantothenic Acid -> <b>pnto__R</b>
#         * 6. Pyridoxine -> <b>pydxn</b>
#         * 7. Biotin -> <b>btn</b>
#         * 9. Folic Acid -> ... can't find anything
#         * 12. Cobalamin -> <b>b12</b> (or cbl1)
# * [Tryptone (2.5) (Oxoid)](http://www.oxoid.com/UK/blue/prod_detail/prod_detail.asp?pr=LP0042&c=UK&lang=EN): 
#     * Tryptophan: <b>trp__L</b>
# * NaCL (1.5) (VWR International, Darmstadt, Germany):
#     * NaCl: <b>na1, cl</b>
# * K<sub>2</sub>HPO<sub>4</sub> (1.0) (Merck International, Darmstadt, Germany):
#     * Potassium: <b>k</b>
#     * Phosphate: <b>p1</b>
#     * Hydrogen: <b>h2</b>
# * KH<sub>2</sub>PO<sub>4</sub> (1.0) (Merck International, Darmstadt, Germany):
#     * Potassium: <b>k</b>
#     * Phosphate: <b>p1</b>
#     * Hydrogen: <b>h2</b>
# * Na<sub>2</sub>SO<sub>4</sub> (2.0) (VWR):
#     * Sodium: <b>na1</b>
#     * Sulfate: <b>so4</b>
# * MgSO<sub>4</sub>*7H<sub>2</sub>O (1.0) (Merck):
#     * Magnesium: <b>mg2</b>
#     * Sulfate: <b>so4</b>
# * CaCl<sub>2</sub>*2H<sub>2</sub>O (0.1) (Merck):
#     * Calcium Chloride: <b>ca2, cl</b>
# * NH<sub>4</sub>Cl (1.0) (Merck):
#     * Ammonium Chloride: <b>nh4, cl</b>
# * Cysteine-HCL (0.4) (Merck):
#     * Cysteine: <b>cys__L</b>
#     * HCL: <b>h2, cl</b>
# * NaHCO<sub>3</sub> (0.2) (VWR):
#     * Sodium: <b>na1</b>
#     * Bicarbonate: <b>hco3<b>
# * MnSO<sub>4</sub>*H<sub>2</sub>O (0.05) (VWR):
#     * Manganese: <b>mn2</b>
#     * Sulfate: <b>so4</b>
# * FeSO<sub>4</sub>*7H<sub>2</sub>O (0.005) (Merck):
#     * Iron: <b>fe</b>
#     * Sulfate: <b>so4</b>
# * ZnSO<sub>4</sub>*7H<sub>2</sub>O (0.005) (VWR):
#     * Zinc: <b>zn2</b>
#     * Sulfate: <b>so4</b>
# * Hemin (0.005) (Sigma-Aldrich, Steinheim, Germany):
#     * ... contains iron
# * Menadione (0.005) (S-A):
#     * <b>mndn</b>
# * Resazurin (0.001) (S-A):
#     * Fluoro identifier
# 
# 
# %% [markdown]
# | Component                            | Concentration (g/L) | BiGG Metabolites                                |
# |--------------------------------------|---------------------|-------------------------------------------------|
# | Bacteriological Peptome              | 6.5                 | polypep, n2                                     |
# | Soy Peptome                          | 5.0                 | stys, raffin, sucr, n2                          |
# | Yeast Extract                        | 3.0                 | thm, ribflv, trp\__L, pnto\__R, pydxn, btn, b12 |
# | Tryptone                             | 2.5                 | trp\__L                                         |
# | NaCL                                 | 1.5                 | na1, cl                                         |
# | K<sub>2</sub>PO<sub>4</sub>          | 1.0                 | k, p1, h2                                       |
# | KH<sub>2</sub>PO<sub>4</sub>         | 1.0                 | k, p1, h2                                       |
# | Na<sub>2</sub>SO<sub>4</sub>         | 2.0                 | na1, so4                                        |
# | Mg<sub>2</sub>SO<sub>4</sub>         | 1.0                 | mg2, so4                                        |
# | CaCl<sub>2</sub>*2H<sub>2</sub>O     | 0.1                 | ca2, cl                                         |
# | NH<sub>4</sub>Cl                     | 1.0                 | nh4, cl                                         |
# | Cysteine-HCL                         | 0.4                 | cys__L, h2, cl                                  |
# | NaHCO<sub>3</sub>                    | 0.2                 | na1, hco3                                       |
# | MnSO<sub>4</sub>*H<sub>2</sub>O      | 0.05                | mn2, so4                                        |
# | FeSOMnSO<sub>4</sub>*7H<sub>2</sub>O | 0.005               | fe, so4                                         |
# | ZnSO<sub>4</sub>*7H<sub>2</sub>O     | 0.005               | zn2, so4                                        |
# | Hemin                                | 0.005               |                                                 |
# | Menadione                            | 0.005               | mndn                                            |
# | Resazurin                            | 0.001               |                                                 |

# %%
model.copy()


# %%
medium_mmcb = list(('polypep','n2','stys','raffin','sucr','n2','thm','ribflv','trp__L',
                    'pnto__R','pydxn','btn','b12','trp__L','na1','cl','k','p1','h2','na1',
                    'so4','mg2','ca2','nh4','cys__L','h2','cl','na1','hco3','mn2','zn2','mndnb'))


# %%
medium_mmcb=[mets+"(e)" for mets in medium_mmcb]
model_mmcb=medium(model.copy(), medium_mmcb)
unused = []
for i in range(1,len(model_mmcb.reactions)):
    model_mmcb.reactions[i].upper_bound=1000
    if model_mmcb.reactions[i].lower_bound != 0:
        model_mmcb.reactions[i].lower_bound=-10
    else:    
        model_mmcb.reactions[i].lower_bound=-
model_mmcb.objective = {model_mmcb.reactions.get_by_id('Biomass_FP'): 1}
model_mmcb.summary()

# Find list of all metabolites not in media
# Express in mmol/g/L
# Set as -'ve flux lower bound


# %%
# Post-optimization
fba_mmcb = model_mmcb.optimize()
fluxes_mmcb = deleteNull(fba_mmcb.fluxes)
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes_mmcb),len(fluxes_mmcb)/len(model.reactions)*100))
print(model_mmcb.objective.value)


