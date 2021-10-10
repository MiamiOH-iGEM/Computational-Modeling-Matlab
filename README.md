# Computational-Modeling-Matlab
All analyses were done using Matlab R2020a and Python3.

**iGEM_Scripts**
All matlab scripts and mat files for the glycolaldehyde model

newglycoaldehyde_gurobi.mat: glycolaldehyde model matlab file

FVAGrowthoneatatime_iGEM.mlx

    Used to perform Flux variability analysis on the models incorporating or not incorporating media, growth, and light data
    
GraphingFVA_iGEM.mlx

    Used to graph and visualize flux variability analysis data saved as csv files from FVAGrowthoneatatime_iGEM.mlx
    
Growthoneatatime_iGEM.mlx

    Used to measure growth and flux of particular reactions over time incorporating media, growth, and light data, saving data as csv files
    
GraphingFluxData_iGEM.mlx

    Used to graph and visualize flux and growth data saved as csv files from Growthoneatatime_iGEM.mlx
    
iGEM_Reaction_List.xlsx

    Excel file showing the added metabolites, reactions, and stoichiometric matrices for the creation of each model (for C4 Collab as well)

**C4_Collaboration**

*Note: the metabolites, reactions, and stoichiometric matrices for the creation of each model can be found under iGEM_Scripts with iGEM_Reaction_List.xlsx

main.py: Python file used to edit the models to create the 3 models of interest by adding the appropriate reactions

C4_MiamiOH.mlx: Matlab script used to load, add sink reactions, and run the analysis on the models

The but, mal, and malbut files contain the following for each model

a. json file for Escher visualization of model

b. mat file for matlab use of model

c. text file indicating the applicable reaction changes to the model

**output_model**

Includes all models in both json format (for Escher visualization) and mat format (for matlab)

**Other**

deletedSBP_Model_iJB785_v1.mat attempts to delete SBPase using remove_genes function while deletedSBP_Model_iJB785_v2.mat attempts to delete SBPase using singleGeneDeletion in python. The models can be used in matlab. The edited_Model_iJB785_v1 (or v2).mat follows the same premise.

main.py is the most updated python script used to create all models excluding the models for the C4 collaboration. The attempt to overexpress TALA is also included in this but was not used.

