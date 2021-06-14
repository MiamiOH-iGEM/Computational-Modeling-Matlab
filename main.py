# importing library
import cobra
from cobra import Reaction, Metabolite


def deleteSBP():
    print("deleting SBP reaction ..... ", end = " ")
    sbp_reaction = model.reactions.SBP
    sbp_reaction.knock_out()
    if ((sbp_reaction.lower_bound == 0) and (sbp_reaction.upper_bound == 0)):
        print("DONE")

def addReaction(rID, rName, rDict):
    reaction = Reaction(rID)
    reaction.name = rName
    # reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    # reaction.objective_coefficient = 0.  # this is the default
    reaction.add_metabolites(rDict)
    print('adding: ', reaction.reaction)
    model.add_reaction(reaction)

#   main function
if __name__ == '__main__':
    global model
    model = cobra.io.load_matlab_model('Model_iJB785_noSpace.mat')
    print('original model: ', model.optimize())
    deleteSBP()

    # Glyceraldehyde_3_phosphate    ------------ R1
    bpg = Metabolite(
        'bpg',
        formula='C3H8O10P2',
        name='1,3_bisphosphoglycerate',
        compartment= 'c')
    dictR1 = {model.metabolites.get_by_id("g3p_c"): -1, bpg: 1}
    addReaction('R1', 'Glyceraldehyde_3_phosphate', dictR1)



    ###################### ignore this part ################################################
    # print(model.metabolites.get_by_id('bpg').formula)

    print('edited model: ', model.optimize())

    # # print(sbp_reaction.genes)
    # gene0505 = model.genes.get_by_id("Synpcc7942_0505")



    # import numpy as np
    # for x in model.metabolites:
    #     print('%9s : %s' % (x.id, x.formula))
    # print(len(model.metabolites))