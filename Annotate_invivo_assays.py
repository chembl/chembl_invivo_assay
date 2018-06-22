#!/nfs/panda/chembl/sw/anaconda3/bin/python
# coding: utf-8

# Use pandas python module to view and analyse data
import pandas as pd
import time
import argparse
import sys
import os


"""
Class description.
"""
__doc__ = """
  AIM: This script annotates the in vivo assays with an assay classification, using text pattern matches. 
 INPUT: invivo_assays.csv, Assay_classification.csv
 OUTPUT: annotated_invivo_assays.csv
 """

#################################################
usage = "Usage: %prog [options] -i input_invivo_assays.csv, -assayclass Assay_Classification_table.csv, -o annotated_invivo_assays.csv" + __doc__

# Example command line instruction:
# bsub python Annotate_invivo_assays.py -i invivo_assays.csv --assayclass Assay_classification.csv -o annotated_invivo_assays.csv
##############################################


##########################################################################################################################################################################
# Function to annotate each in vivo assay (ie a row in the 'input_invivo_assays.csv') using the assay classification:
##########################################################################################################################################################################
#Takes in (i) a dataframe of invivo assays that includes columns for 'assay_chemblid' and 'assay_description'
#       (ii) the assay classification table of reference animal models (based on Hock_2008, Hock_2016, Vogel_2013) and other disease / phenotype / toxicity information

def Annotate(assay_df, Hock_df):
     
    # Add empty new columns to the assay_dataframe:
    assay_df2 = assay_df.copy()
    assay_df2['annotated_classl1'] = None
    assay_df2['annotated_classl2'] = None
    assay_df2['annotated_classl3'] = None
    assay_df2['annotated_text_pattern'] = None
    assay_df2['annotated_source'] = None
    assay_df2['annotated_mesh'] = None
    assay_df2['ann_ind_phenotype_cmpd'] = None
    assay_df2['ann_ind_phenotype_cmpd_chemid'] = None
    assay_df2['ann_pos_control_cmpd'] = None
    assay_df2['ann_pos_control_cmpd_chemblid'] = None
    
    # For each row in the assay classification table, check for pattern matching, or not:
    for row in Hock_df.itertuples():
        annot_l1 = row.l1
        annot_l2 = row.l2
        annot_l3 = row.l3
        annot_text_pattern = row.text_pattern_match
        annot_source = row.assay_source
        annot_mesh = str(row.mesh_indication)
        annot_ind_phenotype_cmpd = str(row.ind_phenotype_cmpd)
        annot_ind_phenotype_cmpd_chemblid = str(row.ind_phenotype_cmpd_chemblid)
        annot_pos_control_cmpd = str(row.pos_control_cmpd)
        annot_pos_control_cmpd_chemblid = str(row.pos_control_cmpd_chemblid)
        
        #If there's a text pattern in the assay classification table, then match against the ChEMBL assay description:
        if pd.notnull(annot_text_pattern):
            subset = assay_df2[ assay_df2['assay_description'].str.contains(annot_text_pattern, case=False)]
#             print("\n", annot_assay_pattern, "; Subset is: ",subset[['annotated_assay_source', 'annotated_assay_mesh']])

            #For each row ('r') that is matched: 
            for r in subset.index:
#                 print("r is:",r,"patt is:", assay_df2.loc[r,'annotated_assay_pattern'])
#                 print("induced_ptype_cmpd is:", assay_df2.loc[ r,'ann_ind_phenotype_cmpd'], "poscontrol is:", type(assay_df2.loc[ r,'ann_pos_control_cmpd']) )
                
                # If there is an existing annotation then append the new reference_assay:
                if assay_df2.loc[ r,'annotated_text_pattern'] != None :
                    assay_df2.at[ r,'annotated_classl1'] = "|".join( [assay_df2.loc[ r,'annotated_classl1'], annot_l1 ] )
                    assay_df2.at[ r,'annotated_classl2'] = "|".join( [assay_df2.loc[ r,'annotated_classl2'], annot_l2] )
                    assay_df2.at[ r,'annotated_classl3'] = "|".join( [assay_df2.loc[ r,'annotated_classl3'], annot_l3] )
                    assay_df2.at[ r,'annotated_text_pattern'] = "|".join( [assay_df2.loc[ r,'annotated_text_pattern'], annot_text_pattern] )
                    assay_df2.at[ r,'annotated_source'] = "|".join( [assay_df2.loc[ r,'annotated_source'], annot_source] )                       
                #Else if there is no existing annotation then add one: 
                elif assay_df2.loc[ r,'annotated_text_pattern'] == None : 
                    assay_df2.at[ r, ['annotated_classl1','annotated_classl2','annotated_classl3'
                                      ,'annotated_text_pattern', 'annotated_source' ]] = [annot_l1,annot_l2, annot_l3, annot_text_pattern, annot_source]
 
                # If there's an existing MeSH_annotation, then append the new reference assay mesh term:
                # Note that need to treat MeSH terms separately since may or may not have a MeSH term for each animal model or disease / phenotype in the assay classification table!!
                if assay_df2.loc[ r,'annotated_mesh'] != None:
                    assay_df2.at[ r,'annotated_mesh'] = "|".join( [assay_df2.loc[ r,'annotated_mesh'], annot_mesh] )
                # Else if there's no existing MeSH_annotation, then add the new reference assay mesh term:    
                elif assay_df2.loc[ r, 'annotated_mesh'] == None:
                    assay_df2.at[ r, 'annotated_mesh'] = annot_mesh

                # Not all reference assays have an annotated compound that induces the phenotype: if present include this too:
                # Note that 'None' is present, and cannot JOIN a NoneType, so need to treat separately:
                if assay_df2.loc[ r,'ann_ind_phenotype_cmpd'] != None:
                    assay_df2.at[ r,'ann_ind_phenotype_cmpd'] = "|".join( [assay_df2.loc[ r,'ann_ind_phenotype_cmpd'], annot_ind_phenotype_cmpd ] ) #Needs str type since includes 'None'
                    assay_df2.at[ r,'ann_ind_phenotype_cmpd_chemid'] = "|".join( [assay_df2.loc[ r,'ann_ind_phenotype_cmpd_chemid'], annot_ind_phenotype_cmpd_chemblid ] ) #Needs str type since includes 'None'
                # If no existing term present for induced_phenotype_compound, then can add directly: 
                elif assay_df2.loc[ r, 'ann_ind_phenotype_cmpd'] == None: 
                    assay_df2.at[ r, 'ann_ind_phenotype_cmpd'] = annot_ind_phenotype_cmpd
                    assay_df2.at[ r, 'ann_ind_phenotype_cmpd_chemid'] = annot_ind_phenotype_cmpd_chemblid
                    
                # Not all reference assays have an annotated positive control compound : if present include this too:
                # Note that 'None' is present, and cannot JOIN a NoneType, so need to treat separately
                if assay_df2.loc[ r,'ann_pos_control_cmpd'] != None:
                    assay_df2.at[ r,'ann_pos_control_cmpd'] = "|".join( [assay_df2.loc[ r,'ann_pos_control_cmpd'], annot_pos_control_cmpd ] )
                    assay_df2.at[ r,'ann_pos_control_cmpd_chemblid'] = "|".join( [assay_df2.loc[ r,'ann_pos_control_cmpd_chemblid'], annot_pos_control_cmpd_chemblid ] )
                # If no existing term present for positive_control_compound, then can add directly: 
                if assay_df2.loc[ r, 'ann_pos_control_cmpd'] == None: 
                    assay_df2.at[ r, 'ann_pos_control_cmpd'] = annot_pos_control_cmpd
                    assay_df2.at[ r, 'ann_pos_control_cmpd_chemblid'] = annot_pos_control_cmpd_chemblid
                elif assay_df2.loc[ r, 'ann_pos_control_cmpd'] == None: 
                    assay_df2.at[ r, 'ann_pos_control_cmpd'] = annot_pos_control_cmpd
                    assay_df2.at[ r, 'ann_pos_control_cmpd_chemblid'] = annot_pos_control_cmpd_chemblid

                    
    return assay_df2
##########################################################################################################################################################################



#################################################
#################################################

#Main section of code:
if __name__ == '__main__':
    
    #Set correct path: 
    path = '/path_to_folder_containing_files/'
#     sys.stdout = open(path+"log_annotate_invivo_assays.txt", 'w') #Redirect std_out to file in output folder!!
    
    #Define parser
    parser = argparse.ArgumentParser(description=__doc__) 
    parser.add_argument("-i", "--invivo_input", dest="invivo_input", default=None, help="input_invivo_assays.csv")
    parser.add_argument("-assayclass", "--assayclass_input", dest="assayclass_input", default=None, help="Assay_Classification_table.csv")
    parser.add_argument("-o", "--output", dest="output", default=None, help="annotated_invivo_assays.csv")
    args = parser.parse_args()  # define arguments that have been parsed        

    print("Run is: " + os.path.basename(sys.argv[0]))
    print("Program starting")
    
    #Read in files for unannotated_invivo_assays & assay_classification_table: 
    invivo_infile = pd.read_csv(path + str(args.invivo_input), skiprows=0 ) 
    assay_class_infile = pd.read_csv(path + str(args.assayclass_input), skiprows=0 ) 
    #Set name of outfile:
    outfile = str(args.output)
    
    ##########################################################################################
    #Now carry out the Annotate function and apply to EACH row of the invivo assays:
    ##########################################################################################

    #Check calculation time:
    start_time = time.time() #Record run time for calculation
    Number_of_assays = len(invivo_infile) # Total number of assays

    df_annotated_invivo_assays = Annotate(invivo_infile, assay_class_infile )

    print(Number_of_assays,"assays in %s minutes" % ( round((time.time() - start_time)/60,2) ))
    ##########################################################################################
    
    ##########################
    # Save table to file:
    df_annotated_invivo_assays.to_csv( path + "/" + outfile, index=False)
    ##########################
    
    print("Program ending")


########################################################################################################################








