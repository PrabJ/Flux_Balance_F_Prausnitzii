"""                    TD PART II :
#################################################################
    STUDY OF THE STRUCTURE OF A GENOME SCALE METABOLIC MODEL
From the given function and the information you got about the cobra language,
try to answer the following question.
files: practical_metaboli_model.pdf
model : study_model.xml
model2 : ####you can chose####
For a better understanding, you can also open the file functions.py to see how some function are built.
During this tutorial you will use object oriented language.
"""
#0 installation : from cmd window 'pip install cobra'

#1- Import needed module
import io
import cobra
import functions
from functions import *

#2- Load the study model
model = io.read_sbml_model('study_model.xml')

#3-	How many metabolites ?
print('Number of metabolites : ',len(model.metabolites))

#4-	How many reactions ?
print('Number of reactions : ',len(model.reactions))

#5-	How many genes 
print('Number of genes : ',len(model.genes))
      
#6-	Which compartments ?
print('Model compartments : ',model.compartments)
      
#7-	Print the stoichiometric matrix
matrix = stoicMatrix(model) #what is this function
print('\nStoichiometric matrix : \n',matrix)
help(stoicMatrix)
#8-	From the stoechiometric matrix , draw the model.
#cf pptx
      
#9-	Identify the external reaction, the transport reaction,  manually and via a command line.
"""
RXNS	TYPE	Cmd line 
Trp_s1	transport	
Trp_s2	Transport	
Trp_p1	Transport	
Ex_s2_e	Exchange	
Ex_s1_e	Exchange	
Ex_p1_e	Exchange	
Rxn_1	Metabolic	
Rxn_2	Metabolic	
Rxn_3	Metabolic	
biomass	biomass	"""
print('Demands : ',model.demands)
print('Exchanges',model.exchanges)
print('Sinks',model.sinks)
print('\nReaction types : \n',reactionsType(model))

#10-	Select a reaction. Print its information. It is a reversible reaction ? Print the bound of the reaction. 
rxn = model.reactions.Rxn_1
print('\nId : %s\nname :%s\nReactionformula : %s\nGenes :%s' %(rxn.id,rxn.name,rxn.reaction,''.join((i.id for i in rxn.genes))))
print('Reversibility : %s' %(rxn.reversibility))
print(rxn.lower_bound, "< Reaction flux value < ",rxn.upper_bound)

#11-	Print a metabolite formula, name Print reactions which involve à specific metabolite (and if possible the reaction associated).
met = model.metabolites.s1_c 
print('\nMetabolite : %s\nName :  %s\nFormula : %s\nCharge : %i\nCompartment : %s' % (met.id,met.name,met.formula,met.charge,model.compartments[met.compartment]))
print('Involved in reactions : \n%s' %('\n'.join(getRxnsFromMet(model,met))))

#12-	What is the objective function ? Can you show the objective coefficient for each reactions?
print('\nObjective function : ', model.objective)
print('\nObjective coefficient:')
for rxn in model.reactions:
    print('%s ->  %i' %(rxn.id,rxn.objective_coefficient))

#####B  Answer the same question with a model of your taste on Bigg Models
import io
import cobra
from functions import *
BiGG_model = io.read_sbml_model('iAM_Pf480.xml')

#13- Donload a model from bigg models database
print("Number of metabolites" , len(BiGG_model.metabolites))
#14- Import the model:

#15- How many reactions?

#16- how many metabolites?

#17- Wich compartments?

#18- Stoichiometric Matrix (show limited at 500/500):

#19- How many exchange reactions ?

#20- How many sink reactions?

#21- How many demand reactions?

#22- What is the objective function ? Maximise or minimize ?


