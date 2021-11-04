import os
import pandas as pd
from script_utils import show_output

def load_primer_df(primer_excel_file):
    '''
    load all files from the excel_file and store in df
    '''
    
    # check for file
    if not os.path.isfile(primer_excel_file):
        show_output(f"Could not find primer database file {primer_excel_file}")
        return
    
    primer_dfs = []    
    for sheet in pd.ExcelFile(primer_excel_file).sheet_names:
        df = pd.read_excel(primer_excel_file, sheet_name=sheet).rename({'status':'Status'}, axis=1)
        df.loc[:, "excel_sheet"] = sheet
        primer_dfs.append(df)
    primer_df = pd.concat(primer_dfs)
    return primer_df


def check_primerDB(mut_df, primer_DB="", padding=25):
    '''
    load a mutation file and find primers in database
        mutations: a pandas dataframe with columns Chr Start End to the minimum
        --> appends a column hits 
    '''
    
    def check_primerDB_row(mut_row, primer_DB_df=pd.DataFrame(), padding=25):
        '''
        checks one mutation entry for existence in primer_df
        '''
                   
        # mark non_existing primer_DB_df with -1 hits
        if primer_DB_df.empty:
            mut_row['hits'] = -1
            return mut_row
        # get the info from the mutation
        start = mut_row['Start'] + padding
        end = mut_row['End'] - padding
        chrom = mut_row['Chr']
        # find overlapping primer entries
        result = primer_DB_df.query('Chr == @chrom and Start < @start and End > @end')
        if (l := len(result.index)):
            # transfer the primer coords to InsertRange
            result.loc[:, 'InsertRange'] = result['Chr'] + ":" + result['Start'].astype(str) + "-" + result['End'].astype(str)
            result.loc[:, 'offsetL'] = mut_row['Start'] - result['Start']
            result.loc[:, 'InsertSize'] = result['End'] - result['Start']
            result.loc[:, 'AmpliconSize'] = result['InsertSize'] + result['fwdPrimer'].str.len() + result['revPrimer'].str.len()
            result.loc[:, 'AmpliconRange'] = result['Chr'] + ":" + (result['Start']- result['fwdPrimer'].str.len()).astype(str) + "-" + (result['End'] + result['revPrimer'].str.len()).astype(str)
            result.loc[:, 'offsetR'] = result['End'] - mut_row['End']
            # 
            for col in ['Chr', 'Start', 'End']:
                result.loc[:,col] = mut_row[col]
            DB_results.append(result.loc[:, ['Chr', 'Start', 'End', 'Plate', 'Pos', 'fwdPrimer', 'revPrimer',
               'Status', 'Temp', 'AmpliconRange', 'AmpliconSize', 'InsertRange',
               'InsertSize', 'offsetL', 'offsetR']])
        mut_row['DBhits'] = l
        return mut_row
    
    if not primer_DB:
        show_output("No primer database file has been provided", color="warning")
        return
    try:
        primer_df = load_primer_df(primer_DB)
        show_output(f"Primer database loaded from {primer_DB} - {len(primer_df.index)} entries found.")
    except:
        show_output("Could not load primer database", color="warning")
    
    
    for col in ['Chr', 'Start', 'End']:
        if not col in mut_df.columns:
            show_output(f"Required column {col} not found in mutation dataframe! Exiting.", color="warning")
            return mut_df
    
    # init empty DB results df
    DB_results = []
    hits_df = mut_df.apply(check_primerDB_row, primer_DB_df=primer_df, padding=-10, axis=1)
    if (hit_count := len(hits_df.query("DBhits > 0"))):
        show_output(f"Found {hit_count} primers in primer database", color="success")
    else:
        show_output("No hits found in primer database", color="normal")
    DB_results_df = pd.concat(DB_results).drop_duplicates(['Chr', 'Start', 'End', 'Plate', 'Pos'])
    return hits_df, DB_results_df