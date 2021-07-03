# importing library
import textwrap

import cobra
from cobra import Reaction, Metabolite
# import computeResult
from cobra.flux_analysis import single_reaction_deletion, flux_variability_analysis

import saveMat
import sys



def deleteSBP():
    print("deleting SBP reaction: ", end= " ")
    print(len(model.reactions), end=" ----> ")
    sbp_reaction = model.reactions.get_by_id("SBP")
    # single_reaction_deletion(model, [sbp_reaction])
    # sbp_reaction.knock_out()
    model.remove_reactions([sbp_reaction])
    # if ((sbp_reaction.lower_bound == 0) and (sbp_reaction.upper_bound == 0)):
    print(len(model.reactions), end=" ..... ")
    print("DONE")

def addReaction(rID, rName, rDict):
    reaction = Reaction(rID)
    reaction.name = rName
    # reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    # reaction.objective_coefficient = 0.  # this is the default
    reaction.add_metabolites(rDict)
    model.add_reaction(reaction)
    # printing result
    printReaction = "adding " + reaction.name + ": " + reaction.reaction
    print('{:<50}{}'.format(printReaction, "========== FBA: "), end="")
    ###### STEP 10 #####
    with model:
        model.objective = reaction.name
        print('{:<20}'.format(model.slim_optimize()))

def computeResult():
    ###### STEP 16  -- Exchange Reaction #####
    print('\n\n ---------------- STEP 16 Exchange Reaction -------------')
    print("Number of Exchange Reaction: ", len(model.exchanges))

    ###### STEP 17  -- Exchange Reaction #####
    print('\n\n ---------------- STEP 17 Consisten/Inconsistent Reaction -------------')
    consistent_model = cobra.flux_analysis.fastcc(model)
    print('Number of consistent reactions: ', len(consistent_model.reactions))

#   main function
if __name__ == '__main__':
    # Write result to a file
    orig_stdout = sys.stdout
    f = open('out.txt', 'w')
    sys.stdout = f


    global model
    model = cobra.io.load_matlab_model('Model_iJB785_noSpace.mat')
    print('original model: ', model.optimize())
    deleteSBP()
    # saveMat.new_save_matlab_model(model, 'deletedSBP_Model_iJB785.mat')

    #### rTCA Pathway
    print('\n\n ---------------- rTCA Pathway -------------')
    # R1
    dpg13_c = Metabolite(
        '13dpg_c',
        formula='C3H8O10P2',
        name='1,3-phospho-d-glyceroyl-phosphate',
        compartment='c')
    dictR1 = {model.metabolites.get_by_id("g3p_c"): -1, dpg13_c: 1}
    addReaction('R1', 'R1', dictR1)

    # R2
    pg3_c = Metabolite(
        '3pg_c',
        formula='C3H7O7P',
        name='3_phosphoglycerate',
        compartment='c')
    dictR2 = {dpg13_c: -1, pg3_c: 1}
    addReaction('R2', 'R2', dictR2)

    # R3
    pg2_c = Metabolite(
        '2pg_c',
        formula='C3H7O7P',
        name='2_phosphoglycerate',
        compartment='c')
    dictR3 = {pg3_c: -1, pg2_c: 1}
    addReaction('R3', 'R3', dictR3)

    # R4
    dictR4 = {pg2_c: -1, model.metabolites.get_by_id("pep_c"): 1}
    addReaction('R4', 'R4', dictR4)

    # R5
    dictR5 = {model.metabolites.get_by_id("oaa_c"): -1,
              model.metabolites.get_by_id("glu__L_c"): -1,
              model.metabolites.get_by_id("akg_c"): 1,
              model.metabolites.get_by_id("asp__L_c"): 1}
    addReaction('R5', 'R5', dictR5)

    # R6
    dictR6 = {model.metabolites.get_by_id("akg_c"): -1,
              model.metabolites.get_by_id("icit_c"): 1}
    addReaction('R6', 'R6', dictR6)

    # R7
    dictR7 = {model.metabolites.get_by_id("icit_c"): -1,
              model.metabolites.get_by_id("cit_c"): 1}
    addReaction('R7', 'R7', dictR7)

    # R8
    dictR8 = {model.metabolites.get_by_id("cit_c"): -1,
              model.metabolites.get_by_id("oaa_c"): 1,
              model.metabolites.get_by_id("accoa_c"): 1}
    addReaction('R8', 'R8', dictR8)

    # R9
    dictR9 = {model.metabolites.get_by_id("accoa_c"): -1,
              model.metabolites.get_by_id("pyr_c"): 1}
    addReaction('R9', 'R9', dictR9)

    # R10
    ddgp_c = Metabolite(
        'ddgp_c',
        formula='C6H8O9P',
        name='2_dehydro_3_deoxy_D_gluconate_6_phosphate',
        compartment='c')
    dictR10 = {model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("pyr_c"): -1,
              ddgp_c: 1}
    addReaction('R10', 'R10', dictR10)

    # R11
    pgc6_c = Metabolite(
        '6pgc_c',
        formula='C6H13O10P',
        name='6_phosphogluconate',
        compartment='c')
    dictR11 = {ddgp_c: -1,
               pgc6_c: 1}
    addReaction('R11', 'R11', dictR11)

    # R12
    glcn_c = Metabolite(
        'glcn_c',
        formula='C6H12O7',
        name='gluconate',
        compartment='c')
    dictR12 = {pgc6_c: -1,
               glcn_c: 1}
    addReaction('R12', 'R12', dictR12)

    # R13
    dictR13 = {glcn_c: -1,
               model.metabolites.get_by_id("ru5p__D_c"): 1}
    addReaction('R13', 'R13', dictR13)

    # R14
    dictR14 = {model.metabolites.get_by_id("ru5p__D_c"): -1,
               model.metabolites.get_by_id("rb15bp_c"): 1}
    addReaction('R14', 'R14', dictR14)


    ############## Revalidatiing after edition ##############
    computeResult()

    # Saving new model
    # saveMat.new_save_matlab_model(model, 'edited_Model_iJB785.mat')


    ###################### ignore this part ################################################
    # print(model.metabolites.get_by_id('bpg').formula)
    print('edited model: ', model.optimize())
    # print("medium", model.medium)
    # # print(sbp_reaction.genes)
    # gene0505 = model.genes.get_by_id("Synpcc7942_0505")


    ###################### Close the file ##################################################
    sys.stdout = orig_stdout
    f.close()




    # import numpy as np
    # for x in model.metabolites:
    #     print('%9s : %s' % (x.id, x.formula))
    # print(len(model.metabolites))