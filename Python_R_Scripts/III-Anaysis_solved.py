"""                    TD PART III :
#################################################################
    ANALYSIS AND HYPOTHESIS WHITH A GENOME SCALE METABOLIC MODEL
From the given functions, the information you got about the cobra language,
and the result of previous TD part, try to answer the following questions.
files: practical_metaboli_model.pdf
model : e_coli_core.xml
        other to download
For a better understanding, you can also open the file functions.py to see how some function are built.
During this tutorial you will use object oriented language.
"""
######A
# 1- Load module:
from functions import *
import cobra
#2- Import e_coli_core model:
cobra.io.save_json_model()
model = cobra.io.read_sbml_model('f_prausnitzii.xml')

#3 - Print information about the model:
info(model)

#4- Lauch the flux balance analysis without changing any parameter:
fba = model.optimize()
#5- Print the objective value. Was it what you expected? 
print(fba.objective_value)

#6- Print fluxes how many reaction are active in this condition?
fluxes = deleteNull(fba.fluxes)
print(fluxes)
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes),len(fluxes)/len(model.reactions)*100))

#7- From the fluxes, define which metabolites are consume and produce by the model?
#Verify your result with the summary function:
"""
PRODUCE     CONSUME
02          co2
glc         h
nh4         h2o
pi
"""
print(model.summary(),'\n')
#8- print the summary of one of this metabolite to see how it's produced and consumed
print(model.metabolites.nh4_c.summary(),'\n')

#9- Verify if the model can grow without oxygene:
## Create a new model to keep initial state

model2 = model.copy()
## Set upper bound and lower bound :
model2.reactions.get_by_id('EX_o2_e').lower_bound = 0

## Launch fba
fba2 = model2.optimize()
fluxes2 = deleteNull(fba2.fluxes)
print(model2.summary())
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes2),len(fluxes2)/len(model.reactions)*100))
## Conclude:
"""YES"""

#10- Verify if the model can use other carbon sources: lactate or formate.
## create a copy of the model
model3 = model.copy()
print(model3.exchanges)
## Close the model
model3= closeModel(model3)
## Change bound to allow entry of the target metabolite.
model3.reactions.get_by_id('EX_lac__D_e').lower_bound =-10
for rxn in model3.exchanges:
    print(rxn.lower_bound,'<',rxn.id,'<',rxn.upper_bound)
fba3 = model3.optimize()
fluxes3 = deleteNull(fba3.fluxes)
print(model3.summary())
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes3),len(fluxes3)/len(model.reactions)*100))


#11- It is possible to produce  acetate and to grow? change the objective function (growth min 0.20):
model4=model.copy()
model4.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM').lower_bound= 0.20
model4.reactions.get_by_id('EX_ac_e').lower_bound= 0.1
model4.objective = {model4.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM'): 1, model4.reactions.get_by_id('EX_ac_e') : 1}
util.solver.linear_reaction_coefficients(model4)
print(model4.objective)
print(model4.summary())
print('\nActive reaction : %i which represent %f%% of the model' %(len(fluxes3),len(fluxes3)/len(model.reactions)*100))
#change the coefficient for ac form 1 to 0.5 for example. What do you observe?


#######B
#Lets look all of this in an entire model;
"""  In biggest model is more difficult to interpret something. Lets see if we can answer to some
question using different functions.
"""

"""  I - You know that you will receive a Saccharomyces cerevisiae strain that you have to grow.
You have 3 different media. Which one will you use to be sure that the strain will grow?"""

# Download one model from bigg database 
model = io.read_sbml_model('iND750.xml')
# Verify that the model can grow without changing any parameter:

print(model.summary())
fba = model.optimize()

# Verify how the model can grow in each medium and conclude:
Medium1 = ['ac_e','ade_e','ala__L_e','arg__L_e','asn__L_e','asp__L_e','ca2_e',
           'cit_e','cys_L_e','man_e','gln__L_e','glu__L_e','h2o_e','his__L_e',
           'ile__L_e','mg2_e','nad_e','nh4_e','nac_e','pi_e','pro__L_e','ser__L_e',
           'so4_e','try__L_e','tyr__L_e','trp__L_e','h2_e','o2_e' ] #c = mannose, - uracile pour medium

model1= medium(model,Medium1)
print('\nMedium1')
print(model1.summary())
Medium2 = ['ac_e','ade_e','ala__L_e','arg__L_e','asn__L_e','asp__L_e','ca2_e',
           'cit_e','cu2_e','cys_L_e','fe2_e','fe3_e','fru_e','glc__D_e',
           'gln__L_e','glu__L_e','h2o_e','his__L_e','ile__L_e','mg2_e','nad_e',
           'nh4_e','nac_e','no3_e','pi_e','pro__L_e','ser__L_e','so4_e','try__L_e',
           'tyr__L_e','ura_e','val__L_e','trp__L_e','h2_e','o2_e' ]# 2 carbon source less mets

model2 = medium(model,Medium2)
print('\nMedium2')
print(model2.summary())
Medium3 = ['4abz_e','ac_e','ace_e','ade_e','ala__L_e','arg__L_e','ascb__L_e','asn__L_e',
           'asp__L_e','btn_e','cbl1_e','cit_e','cl_e','cobalt2_e',
           'cu2_e','cys_L_e','fe2_e','fe3_e','fol_e','glc__D_e','glu__L_e','gly_e','gthrd_e',
           'gua_e','h2o_e','hco3_e','his__L_e','i_e','ile__L_e','inost_e','k_e','leu_L_e','met__L_e','co2_e','no2_e',
           'mg2_e','mn2_e','mndn_e','mobd_e','na1_e','nad_e','nh4_e','nac_e','ni2_e',
           'no3_e','pi_e','pnto__R_e','pro__L_e','pydam_e','pydxn_e','rezrn_e','ribflv_e',
           'ser__L_e','so4_e','thci_e','thm_e','try__L_e','tyr__L_e','ura_e','val__L_e',
           'xan_e','zn2_e','h_e','trp__L_e', 'mobd_e','slnt_e', 'tungs_e','lipoate_e',
           'pydx_e','o2_e' ]# - fructose - gl rich with one carbone source
model3 = medium(model,Medium3)
print('\nMedium3')
help(cobra.medium)
print(model3.summary())

""" II- You work for an agribusiness company specialized in yoghurt production.
In, purpose to improve the taste of some product, you want to produce the molecule RR 2 3 Butanediol which gave the butter/ creamy taste.
Your company whant a biological way to produce this molecule instead of chemical.
you have 4 species in your bacteria collection:
    - Saccharomyces cerevisiae
    - Lactococcus lactis
    - Klebsiella pneumoniae
Can you tell which species is the most likely to produce this molecule? """

# Find molecule id on Bigg database / its exchange reaction associated
met = 'btd_rr_e'
rxn = 'EX_btd_RR_e'

# Load the 3 models

model1 = io.read_sbml_model('iMM904.xml')
model2 = io.read_sbml_model('iNF517.xml')
model3 = io.read_sbml_model('iYL1228.xml')

# Run following function for 'lactococcus lactis':
model2 = updateLactococcus(model2)

# Verify that models can grow:
print('\nSaccharomyces')
print(model1.summary())
print('\nLactococcus')
print(model2.summary())
print('\nKlebsiella')
print(model3.summary())

#for each model find biomass function position or name
type1 = reactionsType(model1)
biomPos1 = type1.index('biomass')

type2 = reactionsType(model2)
biomPos2 = type2.index('biomass')

type3 = reactionsType(model3)
biomPos3 = type3.index('biomass')
# Choice 1 :
#- Update the biomasse bound (>0)
#- Change the objective function
#- Run FBA
print('\nChoice1')
model1.reactions[biomPos1].lower_bound =0.1
model2.reactions[biomPos2].lower_bound =0.1
model3.reactions[biomPos3].lower_bound =0.1

model1.objective = {model1.reactions.get_by_id('EX_btd_RR_e'): 1}
model2.objective = {model2.reactions.get_by_id('EX_btd_RR_e'): 1}
model3.objective = {model3.reactions.get_by_id('EX_btd_RR_e'): 1}

print('\nSaccharomyces')
print(model1.summary())
print('\nLactococcus')
print(model2.summary())
print('\nKlebsiella')
print(model3.summary())

#Choice 2 :
#- change the objective function with both reactions
#- Run FBA
print('\nChoice2')
model1.objective = {model1.reactions.get_by_id('EX_btd_RR_e'): 1, model1.reactions[biomPos1] : 1}
model2.objective = {model2.reactions.get_by_id('EX_btd_RR_e'): 1, model2.reactions[biomPos2] : 1}
model3.objective = {model3.reactions.get_by_id('EX_btd_RR_e'): 1, model3.reactions[biomPos3] : 1}

print('\nSaccharomyces')
print(model1.summary())
print('\nLactococcus')
print(model2.summary())
print('\nKlebsiella')
print(model3.summary())

# What can you conclude?
"""Lactobacilus lactis is the best producer of this molecule with th highet flux in the exchange
reaction."""

# Did you check if models could spontaneously produce the molecule?
""" Lactococcus lactis is producing the molecule when we launch the FBA for growth.
we think that its more likly to produce the molecule"""


