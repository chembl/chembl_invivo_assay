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
  AIM: This script annotates the in vivo assays with Hock reference animal models, using text pattern matches. 
 INPUT: invivo_assays.csv, Hock_table.csv
 OUTPUT: annotated_invivo_assays.csv
 """

#################################################
usage = "Usage: %prog [options] -i input_invivo_assays.csv, -hock Hock_table.csv, -o annotated_invivo_assays.csv" + __doc__
# Example command line instruction:
# bsub python Annotate_invivo_assays_commandline.py -i invivo_assays.csv -hock Hock_reference_assays.csv -o annotated_invivo_assays.csv
##############################################


##########################################################################################################################################################################
# Function to annotate each in vivo assay (ie a row in the 'input_invivo_assays.csv') using the Hock reference animal models:
##########################################################################################################################################################################
#Takes in (i) a dataframe of invivo assays that includes columns for 'assay_chemblid' and 'assay_description'
#       (ii) the Hock_Vogel table of reference animal models

def Annotate(assay_df, Hock_df):
     
    # Add empty new columns to the assay_dataframe:
    assay_df2 = assay_df.copy()
    assay_df2['annotated_chaptercode'] = None
    assay_df2['annotated_chapter'] = None
    assay_df2['annotated_subchapter'] = None
    assay_df2['annotated_assay_refname'] = None
    assay_df2['annotated_assay_pattern'] = None
    assay_df2['annotated_assay_source'] = None
    assay_df2['annotated_assay_mesh'] = None
    assay_df2['ann_ind_phenotype_cmpd'] = None
    assay_df2['ann_ind_phenotype_cmpd_chemid'] = None
    assay_df2['ann_pos_control_cmpd'] = None
    assay_df2['ann_pos_control_cmpd_chemblid'] = None
    
    # For a row in the Hock_df animal models table, check for text pattern matching, or not:
    for row in Hock_df.itertuples():
        annot_assay_chptcode = row.chapter_code
        annot_assay_chapter = row.chapter
        annot_assay_subchapter = row.subchapter
        annot_assay_refname = row.assay_refname
        annot_assay_pattern = row.text_pattern_match
        annot_assay_source = row.assay_source
        annot_assay_mesh = str(row.mesh_indication)
        annot_ind_phenotype_cmpd = str(row.ind_phenotype_cmpd)
        annot_ind_phenotype_cmpd_chemblid = str(row.ind_phenotype_cmpd_chemblid)
        annot_pos_control_cmpd = str(row.pos_control_cmpd)
        annot_pos_control_cmpd_chemblid = str(row.pos_control_cmpd_chemblid)
        
        #If there's a pattern in the Hock_df, then match against the invivo assay:
        if pd.notnull(annot_assay_pattern):
            subset = assay_df2[ assay_df2['assay_description'].str.contains(annot_assay_pattern, case=False)]
#             print("\n", annot_assay_pattern, "; Subset is: ",subset[['annotated_assay_source', 'annotated_assay_mesh']])

            #For each row ('r') that is matched: 
            for r in subset.index:
#                 print("r is:",r,"patt is:", assay_df2.loc[r,'annotated_assay_pattern'])
#                 print("induced_ptype_cmpd is:", assay_df2.loc[ r,'ann_ind_phenotype_cmpd'], "poscontrol is:", type(assay_df2.loc[ r,'ann_pos_control_cmpd']) )
                
                # If there is an existing annotation then append the new reference_assay:
                if assay_df2.loc[ r,'annotated_assay_pattern'] != None :
                    assay_df2.at[ r,'annotated_chaptercode'] = "|".join( [assay_df2.loc[ r,'annotated_chaptercode'], annot_assay_chptcode ] )
                    assay_df2.at[ r,'annotated_chapter'] = "|".join( [assay_df2.loc[ r,'annotated_chapter'], annot_assay_chapter ] )
                    assay_df2.at[ r,'annotated_subchapter'] = "|".join( [assay_df2.loc[ r,'annotated_subchapter'], annot_assay_subchapter] )
                    assay_df2.at[ r,'annotated_assay_refname'] = "|".join( [assay_df2.loc[ r,'annotated_assay_refname'], annot_assay_refname] )
                    assay_df2.at[ r,'annotated_assay_pattern'] = "|".join( [assay_df2.loc[ r,'annotated_assay_pattern'], annot_assay_pattern] )
                    assay_df2.at[ r,'annotated_assay_source'] = "|".join( [assay_df2.loc[ r,'annotated_assay_source'], annot_assay_source] )                       
                #Else if there is no existing annotation then add one: 
                elif assay_df2.loc[ r,'annotated_assay_pattern'] == None : 
                    assay_df2.at[ r, ['annotated_chaptercode', 'annotated_chapter','annotated_subchapter'
                                      ,'annotated_assay_refname','annotated_assay_pattern', 'annotated_assay_source' ]] = [annot_assay_chptcode, annot_assay_chapter,annot_assay_subchapter
                                                                                                                        , annot_assay_refname,annot_assay_pattern, annot_assay_source]
 
                # If there's an existing MeSH_annotation, then append the new reference assay mesh term:
                # Note that need to treat MeSH terms separately since may or may not have a MeSH term for each Hock_reference_assay!!
                if assay_df2.loc[ r,'annotated_assay_mesh'] != None:
                    assay_df2.at[ r,'annotated_assay_mesh'] = "|".join( [assay_df2.loc[ r,'annotated_assay_mesh'], annot_assay_mesh] )
                # Else if there's no existing MeSH_annotation, then add the new reference assay mesh term:    
                elif assay_df2.loc[ r, 'annotated_assay_mesh'] == None:
                    assay_df2.at[ r, 'annotated_assay_mesh'] = annot_assay_mesh

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
    path = '/nfs/panda/chembl/shared/hecatos/invivo_annotation/'
#     sys.stdout = open(path+"log_annotate_invivo_assays.txt", 'w') #Redirect std_out to file in output folder!!
    
    #Define parser
    parser = argparse.ArgumentParser(description=__doc__) 
    parser.add_argument("-i", "--invivo_input", dest="invivo_input", default=None, help="input_invivo_assays.csv")
    parser.add_argument("-hock", "--hock_input", dest="hock_input", default=None, help="Hock_table.csv")
    parser.add_argument("-o", "--output", dest="output", default=None, help="annotated_invivo_assays.csv")
    args = parser.parse_args()  # define arguments that have been parsed        

    print("Run is: " + os.path.basename(sys.argv[0]))
    print("Program starting")
    
    #Read in files for unannotated_invivo_assays & Hock_table: 
    invivo_infile = pd.read_csv(path + str(args.invivo_input), skiprows=0 ) 
    hock_infile = pd.read_csv(path + str(args.hock_input), skiprows=0 ) 
    #Set name of outfile:
    outfile = str(args.output)
    
    ##########################################################################################
    #Now carry out the Annotate function and apply to EACH row of the invivo assays:
    ##########################################################################################

    #Check calculation time:
    start_time = time.time() #Record run time for calculation
    Number_of_assays = len(invivo_infile) # Total number of assays

    df_annotated_invivo_assays = Annotate(invivo_infile, hock_infile )

    print(Number_of_assays,"assays in %s minutes" % ( round((time.time() - start_time)/60,2) ))
    ##########################################################################################
    
    ##########################
    # Save table to file:
    df_annotated_invivo_assays.to_csv( path + "/" + outfile, index=False)
    ##########################
    
    print("Program ending")


########################################################################################################################








