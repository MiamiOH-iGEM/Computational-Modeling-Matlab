# importing library
import textwrap

import cobra
from cobra import Reaction, Metabolite
# import computeResult
from cobra.flux_analysis import single_reaction_deletion, flux_variability_analysis, single_gene_deletion
from cobra.manipulation import remove_genes

import saveMat
import sys



def deleteSBP():
    initalGeneNumber = len(model.genes)
    print("\n---------------- Deleting the  SBP reaction ----------------")
    print("deleting SBP reaction: ", end= " ")
    print(len(model.reactions), end=" ----> ")
    sbp_reaction = model.reactions.get_by_id("SBP")
    # single_reaction_deletion(model, [sbp_reaction])
    # sbp_reaction.knock_out()     # if ((sbp_reaction.lower_bound == 0) and (sbp_reaction.upper_bound == 0)):
    model.remove_reactions([sbp_reaction])
    print(len(model.reactions), end=" ..... DONE\n")
    deletedSBP_model = model.optimize()
    print('Growth (after deleting SBP): ', deletedSBP_model.objective_value, '\n')

    print("\n---------------- Deleting the  SBP reaction gene (Synpcc7942_0505) ----------------")
    # single_gene_deletion(model, [model.genes.get_by_id("Synpcc7942_0505")])    # simulating gene deletion
    remove_genes(model, sbp_reaction.genes)                                  # deleting the Synpcc7942_0505 gene
    print("Number of genes: ", initalGeneNumber, "---->", len(model.genes), "..... DONE")
    deletedSBPgene_model = model.optimize()
    print("Growth (after deleting SBP's gene): ", deletedSBPgene_model.objective_value, '\n')

def addReaction(rID, rName, rDict, geneStr = None):
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
    # print('{:<50}{}'.format(printReaction, "========== FBA: "), end="")

    # Adding gene rule
    if (geneStr != None):
        reaction.gene_reaction_rule = geneStr
    print('{:<50}{}{}'.format(printReaction, "--- gene: ", reaction.genes))


    # ###### STEP 10 #####
    # with model:
    #     model.objective = reaction.name
    #     print('{:<20}'.format(model.slim_optimize()))

def computeResult():
    ###### STEP 16  -- Exchange Reaction #####
    print('\n\n---------------- STEP 16 Exchange Reaction -------------')
    print("Number of Exchange Reaction: ", len(model.exchanges))

    ###### STEP 17  -- Exchange Reaction #####
    print('\n\n ---------------- STEP 17 Consisten/Inconsistent Reaction -------------')
    consistent_model = cobra.flux_analysis.fastcc(model)
    print('Number of consistent reactions: ', len(consistent_model.reactions))

def addrTCA():
    print('\n---------------- Adding rTCA Pathway -------------')
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

def addTransaldolaseAlternativePathway():
    print('\n---------------- Adding Transaldolase Alternative Pathway -------------')
    # R1
    dictR1 = {model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("dhap_c"): -1,
              model.metabolites.get_by_id("fdp_c"): 1}
    addReaction('R1', 'R1', dictR1, 'SYNPCC7942_RS07385')

    # R2
    dictR2 = {model.metabolites.get_by_id("fdp_c"): -1,
              model.metabolites.get_by_id("f6p_c"): 1}
    addReaction('R2', 'R2', dictR2, 'SYNPCC7942_RS11675')

    # R3
    dictR3 = {model.metabolites.get_by_id("g3p_c"): 1,
              model.metabolites.get_by_id("f6p_c"): -1,
              model.metabolites.get_by_id("e4p_c"): -1,
              model.metabolites.get_by_id("s7p_c"): 1 }
    addReaction('R3', 'R3', dictR3,'SYNPCC7942_RS10545')

    # R4
    dictR4 = {model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("s7p_c"): -1,
              model.metabolites.get_by_id("r5p_c"): 1,
              model.metabolites.get_by_id("xu5p__D_c"): 1}
    addReaction('R4', 'R4', dictR4,'SYNPCC7942_RS02735')

    # R5
    dictR5 = {model.metabolites.get_by_id("r5p_c"): -1,
              model.metabolites.get_by_id("ru5p__D_c"): 1}
    addReaction('R5', 'R5', dictR5,'SYNPCC7942_RS02970')

    # R6
    dictR6 = {model.metabolites.get_by_id("xu5p__D_c"): -1,
              model.metabolites.get_by_id("ru5p__D_c"): 1}
    addReaction('R6', 'R6', dictR6,'SYNPCC7942_RS03100')

    # R7
    dictR7 = {model.metabolites.get_by_id("ru5p__D_c"): -1,
              model.metabolites.get_by_id("rb15bp_c"): 1}
    addReaction('R7', 'R7', dictR7, 'SYNPCC7942_RS05015')

def addEditedGlycolaldehydeAdditionPathway():
    print('\n---------------- Adding Edited Glycolaldehyde Addition Pathway -------------')
    # R1
    dictR1 = {model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("gcald_c"): -1,
              model.metabolites.get_by_id("ara5p_c"): 1}
    addReaction('R1', 'R1', dictR1)

    # R2
    dictR2 = {model.metabolites.get_by_id("ara5p_c"): -1,
              model.metabolites.get_by_id("ru5p__D_c"): 1}
    addReaction('R2', 'R2', dictR2, '( Synpcc7942_2291 or SYNPCC7942_RS11645 )')


    # R3
    dictR3 = {model.metabolites.get_by_id("ru5p__D_c"): -1,
              model.metabolites.get_by_id("rb15bp_c"): 1}
    addReaction('R3', 'R3', dictR3, '( Synpcc7942_0977 or  SYNPCC7942_RS05015 )')

    # R4
    ru_D_c = Metabolite(
        'ru_D_c',
        formula='C5H10O5',
        name='d-ribulose',
        compartment='c')
    dictR4 = {model.metabolites.get_by_id("ru5p__D_c"): -1,
              ru_D_c: 1}
    addReaction('R4', 'R4', dictR4, '( Synpcc7942_2462 or SYNPCC7942_RS12515 )')

    # R5
    ru1p_c = Metabolite(
        'ru1p_c',
        formula='C5H11O8P',
        name='ribulose-1-phosphate',
        compartment='c')
    dictR5 = {ru_D_c: -1, ru1p_c: 1}
    addReaction('R5', 'R5', dictR5)

    # R6
    dictR6 = {model.metabolites.get_by_id("gcald_c"): 1,
              ru1p_c: -1,
              model.metabolites.get_by_id("dhap_c"): 1}
    addReaction('R6', 'R6', dictR6)

    # R7
    dictR7 = {model.metabolites.get_by_id("g3p_c"): 1,
              model.metabolites.get_by_id("dhap_c"): -1}
    addReaction('R7', 'R7', dictR7, '( Synpcc7942_1261 or SYNPCC7942_RS06460 )')

def addNewEditedGlycolaldehydeAdditionPathway():
    print('\n---------------- Adding Edited Glycolaldehyde Addition Pathway -------------')
    # R1
    dictR1 = {model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("gcald_c"): -1,
              model.metabolites.get_by_id("ara5p_c"): 1}
    addReaction('R1', 'R1', dictR1)

    # R4
    dictR4 = {model.metabolites.get_by_id("g3p_c"): -2,
              model.metabolites.get_by_id("gcald_c"): 3,
              model.metabolites.get_by_id("pi_c"): 2}
    addReaction('R4', 'R4', dictR4)

def addGlycolaldehydeAdditionPathway():
    print('\n---------------- Adding Glycolaldehyde Addition Pathway -------------')
    # R1
    dictR1 = {model.metabolites.get_by_id("2pglyc_c"): -1,
              model.metabolites.get_by_id("pi_c"): 1,
              model.metabolites.get_by_id("glyclt_c"): 1}
    addReaction('R1', 'R1', dictR1, 'Synpcc7942_2613')

    # R2
    glyCoa = Metabolite(
        'glyCoa',
        formula='C23H38N7O18P3S',
        name='d-Glycolyl-CoA',
        compartment='c')
    dictR2 = {model.metabolites.get_by_id("glyclt_c"): -1,
              model.metabolites.get_by_id("coa_c"): -1,
              glyCoa: 1,               ############## check again
              model.metabolites.get_by_id("ppi_c"): 1}
    addReaction('R2', 'R2', dictR2, '6GVS_J')


    # R3
    dictR3 = {model.metabolites.get_by_id("glyald_c"): 1,
              glyCoa: -1 ,            ################### check again
              model.metabolites.get_by_id("coa_c"): 1}
    addReaction('R3', 'R3', dictR3, 'Synpcc7942_0489')

    # R4
    dictR4 = {model.metabolites.get_by_id("glyald_c"): -1,
              model.metabolites.get_by_id("g3p_c"): -1,
              model.metabolites.get_by_id("ara5p_c"): 1}
    addReaction('R4', 'R4', dictR4, 'Synpcc7942_1443')

    # R5
    dictR5 = {model.metabolites.get_by_id("ara5p_c"): -1,
              model.metabolites.get_by_id("ru5p__D_c"): 1}
    addReaction('R5', 'R5', dictR5, 'Synpcc7942_2291')

    # R6
    dictR6 = {model.metabolites.get_by_id("ru5p__D_c"): -1,
              model.metabolites.get_by_id("rb15bp_c"): 1}
    addReaction('R6', 'R6', dictR6, 'Synpcc7942_0977')

def talaOverexpress():
    print("\n---------------- Overexpressing TALA genes/reaction ----------------")
    # ================ Start to overexpress =====================================
    tala_reaction = model.reactions.get_by_id("TALA")
    fva = flux_variability_analysis(model, reaction_list=[tala_reaction], fraction_of_optimum=0.9)
    print("----- Reaction optimized range ------")
    print(fva)

    model.reactions.get_by_id('TALA').lower_bound = -2
    overexpress_model = model.optimize()  # solution is stored at model.solution
    print("Growth (after overexpressing TALA): ", overexpress_model.objective_value)

def fbpOverexpress():
    print("\n---------------- Overexpressing FBP genes/reaction ----------------")
    # ================ Start to overexpress =====================================
    fbp_reaction = model.reactions.get_by_id("FBP")
    fva = flux_variability_analysis(model, reaction_list=[fbp_reaction], fraction_of_optimum=0.9)
    print("----- Reaction optimized range ------")
    print(fva)

    model.reactions.get_by_id('FBP').lower_bound = 1000
    overexpress_model = model.optimize()  # solution is stored at model.solution
    print("Growth (after overexpressing FBP): ", overexpress_model.objective_value)

def addButanolBiosyn():
    print('\n---------------- Adding Butanol Biosynthesis Pathway -------------')
    # # R1
    # dictR1 = {model.metabolites.get_by_id("accoa_c"): -1,
    #           model.metabolites.get_by_id("co2_c"): -1,
    #           model.metabolites.get_by_id("malcoa_c"): 1}
    # addReaction('R1', 'R1', dictR1, 'Synpcc7942_1379')

    # R2
    aAcoa_c = Metabolite(
        'aAcoa_c',
        formula='C25H40N7O18P3S',
        name='AcetoAcetyl-CoA',
        compartment='c')
    dictR2 = {model.metabolites.get_by_id("accoa_c"): -2,
              model.metabolites.get_by_id("coa_c"): 1,
              aAcoa_c: 1}
    addReaction('R2', 'R2', dictR2)

    # R3
    dictR3 = {model.metabolites.get_by_id("accoa_c"): -1,
              model.metabolites.get_by_id("malcoa_c"): -1,
              aAcoa_c: 1,
              model.metabolites.get_by_id("co2_c"): 1}
    addReaction('R3', 'R3', dictR3, 'NphT7')

    # R4
    ohbcoa_c = Metabolite(
        'ohbcoa_c',
        formula='C25H38N7O18P3S',
        name='R-3-OH-Butyril-CoA',
        compartment='c')
    dictR4 = {model.metabolites.get_by_id("h_c"): -1,
              model.metabolites.get_by_id("nadph_c"): -1,
              model.metabolites.get_by_id("nadp_c"): 1,
              ohbcoa_c: 1,
              aAcoa_c: -1}
    addReaction('R4', 'R4', dictR4, 'phaB')

    # R5
    crcoa_c = Metabolite(
        'crcoa_c',
        formula='C25H40N7O17P3S',
        name='Crotonyl-CoA',
        compartment='c')
    dictR5 = {model.metabolites.get_by_id("h2o_c"): 1,
              ohbcoa_c: -1,
              crcoa_c: 1}
    addReaction('R5', 'R5', dictR5, 'Ccr')

    # R6
    butcoa_c = Metabolite(
        'butcoa_c',
        formula='C25H42N7O17P3S',
        name='Butyril_CoA',
        compartment='c')
    dictR6 = {model.metabolites.get_by_id("h_c"): -1,
              model.metabolites.get_by_id("nadh_c"): -1,
              model.metabolites.get_by_id("nad_c"): 1,
              crcoa_c: -1,
              butcoa_c: 1}
    addReaction('R6', 'R6', dictR6, 'PhaJ')

    # R7
    but_c = Metabolite(
        'but_c',
        formula='C4H8O2',
        name='butyrate',
        compartment='c')
    dictR7 = {model.metabolites.get_by_id("coa_c"): 1,
              model.metabolites.get_by_id("h2o_c"): -1,
              butcoa_c: -1,
              but_c: 1}
    addReaction('R7', 'R7', dictR7, 'PduP')

    # R8
    butaldh_c = Metabolite(
        'butaldh_c',
        formula='C4H8O',
        name='Butyraldehyde',
        compartment='c')
    dictR8 = {model.metabolites.get_by_id("atp_c"): -1,
              model.metabolites.get_by_id("amp_c"): 1,
              model.metabolites.get_by_id("pi_c"): 2,
              model.metabolites.get_by_id("nadph_c"): -1,
              model.metabolites.get_by_id("nadp_c"): 1,
              butaldh_c: 1,
              but_c: -1}
    addReaction('R8', 'R8', dictR8, 'Slr1192')

    # R9
    nbut_c = Metabolite(
        'nbut_c',
        formula='C4H10O',
        name='n-Butanol',
        compartment='c')
    dictR9 = {model.metabolites.get_by_id("h_c"): -1,
              model.metabolites.get_by_id("nadh_c"): -1,
              model.metabolites.get_by_id("nad_c"): 1,
              butaldh_c: -1,
              nbut_c: 1}
    addReaction('R9', 'R9', dictR9)

def addMelylCoa():
    print('\n---------------- Adding Malyl-CoA-Glycerate Artificial Cycle + PK Bypass -------------')
    # # R1
    # dictR1 = {model.metabolites.get_by_id("pep_c"): -2,
    #           model.metabolites.get_by_id("oaa_c"): 2,
    #           model.metabolites.get_by_id("hco3_c"): -2}
    # addReaction('R10', 'R10', dictR1, 'ppc')

    # R2
    dictR2 = {model.metabolites.get_by_id("oaa_c"): -1,
              model.metabolites.get_by_id("mal__L_c"): 1,
              model.metabolites.get_by_id("nadh_c"): -1,
              model.metabolites.get_by_id("nad_c"): 1,
              model.metabolites.get_by_id("h_c"): -1}
    addReaction('R20', 'R20', dictR2,'mdh')

    # R3
    malycoa_c = Metabolite(
        'malycoa_c',
        formula='C25H40N7O20P3S',
        name='Malyl-oCA',
        compartment='c')
    dictR3 = {model.metabolites.get_by_id("mal__L_c"): -1,
              malycoa_c: 1,
              model.metabolites.get_by_id("coa_c"): -1,
              model.metabolites.get_by_id("pi_c"): 1,
              model.metabolites.get_by_id("atp_c"): -1,
              model.metabolites.get_by_id("adp_c"): 1}
    addReaction('R30', 'R30', dictR3, 'mtk')

    # R4
    dictR4 = {malycoa_c: -1,
              model.metabolites.get_by_id("glx_c"): 1,
              model.metabolites.get_by_id("accoa_c"): 1}
    addReaction('R40', 'R40', dictR4, 'mcl')

    # R5
    dictR5 = {model.metabolites.get_by_id("glx_c"): -2,
              model.metabolites.get_by_id("2h3oppan_c"): 1,
              model.metabolites.get_by_id("co2_c"): 1}
    addReaction('R50', 'R50', dictR5, 'gcl')

    # # R6
    # dictR6 = {tar_saldh: -1,
    #           model.metabolites.get_by_id("glyc__R_c"): 1}
    # addReaction('R60', 'R60', dictR6, 'garR')

    # # R7
    # dictR7 = {model.metabolites.get_by_id("glyc__R_c"): -1,
    #           model.metabolites.get_by_id("2pg_c"): 1}
    # addReaction('R70', 'R70', dictR7, 'garK')

    # # R8
    # dictR8 = {model.metabolites.get_by_id("pep_c"): 1,
    #           model.metabolites.get_by_id("2pg_c"): -1}
    # addReaction('R80', 'R80', dictR8, 'Eno')

    # R9
    # dictR9 = {model.metabolites.get_by_id("xu5p__D_c"): -1,
    #           model.metabolites.get_by_id("g3p_c"): 1,
    #           model.metabolites.get_by_id("actp_c"): 1}
    # addReaction('R90', 'R90', dictR9, 'PK')

    # R10
    dictR10 = {model.metabolites.get_by_id("pi_c"): 1,
               model.metabolites.get_by_id("accoa_c"): 1,
               model.metabolites.get_by_id("coa_c"): -1,
               model.metabolites.get_by_id("actp_c"): -1}
    addReaction('R100', 'R100', dictR10, 'PtaBs')


#   main function
if __name__ == '__main__':
    # Write result to a file
    orig_stdout = sys.stdout
    f = open('test.txt', 'w')
    sys.stdout = f

    global model
    model = cobra.io.load_matlab_model('Model_iJB785_noSpace.mat')
    model.solver = "gurobi"
    original_model = model.optimize()
    print("Model using: ", type(model.solver))
    print('Growth (original model): ', original_model.objective_value, '\n')
    print('Number of Genes: ', len(model.genes))

    # =============== Saving New Json Origin Model  =========================
    # cobra.io.save_json_model(model, "./outputModel/json/Model_iJB785_noSpace.json")

    # =============== Deleting SBP =========================
    deleteSBP()
    # =============== Saving New Json Model  =========================
    # cobra.io.save_json_model(model, "./outputModel/json/deletedSBP_Model_iJB785_noSpace.json")

    # =============== Saving matlab file of deleting SBP =========================
    # saveMat.new_save_matlab_model(model, 'deletedSBP_Model_iJB785_v2.mat')

    # =============== Adding rTCA Pathway =========================
    # addrTCA()
    # rTCA_model = model.optimize()  # solution is stored at model.solution
    # print("Growth (after adding rTCA pathway): ", rTCA_model.objective_value)

    # =============== Adding Transaldolase ALternative Pathway =========================
    # addTransaldolaseAlternativePathway()
    # transaldolaseAlternative_model = model.optimize()  # solution is stored at model.solution
    # print("Growth (after adding Transaldolase ALternative pathway): ", transaldolaseAlternative_model.objective_value)

    # =============== Adding Edited Glycoaldehyde Addition Pathway =========================
    # addEditedGlycolaldehydeAdditionPathway()
    # editedGlycolaldehydeAddition_model = model.optimize()  # solution is stored at model.solution
    # print("Growth (after adding Edited Glycoaldehyde Addition pathway): ", editedGlycolaldehydeAddition_model.objective_value)

    # =============== Adding New Edited Glycoaldehyde Addition Pathway =========================
    # addNewEditedGlycolaldehydeAdditionPathway()
    # newGlycolaldehydeAddition_model = model.optimize()  # solution is stored at model.solution
    # print("Growth (after adding Edited Glycoaldehyde Addition pathway): ", newGlycolaldehydeAddition_model.objective_value)

    # =============== Adding Glycoaldehyde Addition Pathway =========================
    # addGlycolaldehydeAdditionPathway()
    # glycolaldehydeAddition_model = model.optimize()  # solution is stored at model.solution
    # print("Growth (after adding Glycoaldehyde Addition pathway): ", glycolaldehydeAddition_model.objective_value)

    # =============== Overexpressed genes/reaction =========================
    # talaOverexpress()
    # fbpOverexpress()

    # =============== Revalidating after edition  =========================
    # computeResult()

    # =============== Adding Butanol Biosynthesis Pathway =========================
    addButanolBiosyn()
    butanolBiosyn_model = model.optimize()  # solution is stored at model.solution
    print("Growth (after adding Glycoaldehyde Addition pathway): ", butanolBiosyn_model.objective_value)
    #
    # # =============== Adding Malyl-CoA-Glycerate Pathway =========================
    addMelylCoa()
    malylCoa_model = model.optimize()  # solution is stored at model.solution
    print("Growth (after adding Malyl-CoA-Glycerate pathway): ", malylCoa_model.objective_value)

    # =============== Saving New Matlab Model  =========================
    saveMat.new_save_matlab_model(model, './outputModel/malylCoaAndButan_gurobi.mat')

    # =============== Saving New Json Model  =========================
    cobra.io.save_json_model(model, "./outputModel/json/malylCoaAndButan_gurobi.json")

    # =============== Print Final Growth  =========================
    print('\n\n---------------- Final Model -------------')
    final_model = model.optimize()
    print('Growth (Final) : ', final_model.objective_value)
    print('Number of Genes: ', len(model.genes))

    # =============== Close the file  =========================
    sys.stdout = orig_stdout
    f.close()