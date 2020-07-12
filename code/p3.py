import primer3
import pandas as pd
from primer_utils import mut2insert
from genomics import get_chrom

############# DEFAULT CONFIGS ##################
# for description of config see:
#       manual in https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_manual.htm#PRIMER_MAX_SELF_ANY
PCR_config = {
    'seq_len': 500,
    'mut_pad': 5,
    'prod_size_min': 120,
    'prod_size_max': 220
}

primer3_config = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,  # set 1 to actually use the thermodynamic calculations
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_WT_SELF_END':1,   # use Primer_max_self_end
        'PRIMER_MAX_SELF_END_TH': 30,
        'PRIMER_WT_SELF_END_TH':1, # Primer_max_self_end_th
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 8,
        'PRIMER_PAIR_WT_COMPL_ANY':2,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 30,
        'PRIMER_PAIR_WT_COMPL_ANY_TH':2,
        'PRIMER_MAX_HAIRPIN_TH': 47,
        'PRIMER_WT_HAIRPIN_TH':1
    }

def get_primer_pair(row, chrom_dict={}, config={}):
    '''
    returns the best primer pair for a given location row
    return value is [fwd_seq, fwd_tmp, rev_seq, rev_tmp, prod_size]
    active chromosome sequence is global variable chrom
    '''

    # load sequence
    pos = row['Start']
    half_seq = int(config['seq_len'] / 2)
    seq_start = pos - half_seq
    seq_end = pos + half_seq
    seq = chrom_dict['sequence'][seq_start:seq_end]
    pad = int(config['mut_pad'] / 2)
    half_size = int(config['prod_size_min'] / 2)

    # calculate the target_range as offSet from center (half)
    offSet = half_size - 20 - pad
    target_start = half_seq - offSet
    target = [target_start, offSet * 2]
    setting = {
        'SEQUENCE_ID': 'asdf',
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_TARGET': target
    }
    primers = primer3.bindings.designPrimers(setting, config)

    # return '--' if nothing was found
    if primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
        row['fwd_seq'] = row['fwd_tmp'] = row['rev_seq'] = row['rev_tmp'] = row['prod_size'] = '--'
        return row
    
    ## get chrom coords
    amp_start = seq_start + primers['PRIMER_LEFT_0'][0] + 1
    amp_end = seq_start + primers['PRIMER_RIGHT_0'][0] + 1
    
    row['AmpliconRange'] = f"{chrom_dict['name']}:{amp_start}-{amp_end}"
    
    insert_start = amp_start + primers['PRIMER_LEFT_0'][1]
    insert_end = amp_end - primers['PRIMER_RIGHT_0'][1]
    row['InsertRange'] = f"{chrom_dict['name']}:{insert_start}-{insert_end}"
    row['InsertSize'] = insert_end - insert_start + 1
    insert_seq = chrom_dict['sequence'][insert_start -1:insert_end]
    row['InsertSeq'] = insert_seq
    
    mut_dict = {
        'Chr': row['Chr'],
        'Start': row['Start'],
        'End': row['End'],
        'Ref': row['Ref'],
        'Alt': row['Alt']
    }
    
    insert_dict = {
        'Chr': row['Chr'],
        'Start': insert_start,
        'End': insert_end,
        'seq': insert_seq
    }
    
    
    row['InsertSeq'] = mut2insert(mut=mut_dict, seq=insert_dict, return_none=True)
    row['offsetL'] = row['Start'] - insert_start
    row['offsetR'] = insert_end - row['End']
    
    row['fwdPrimer'] = primers['PRIMER_LEFT_0_SEQUENCE']
    row['revPrimer'] = primers['PRIMER_RIGHT_0_SEQUENCE']
    row['AmpliconSize'] = primers['PRIMER_RIGHT_0'][0] - primers['PRIMER_LEFT_0'][0] + 1
    row['Status'] = 'not established'
    row['Temp'] = f"fwd={int(primers['PRIMER_LEFT_0_TM'] * 10) / 10}|rev={int(primers['PRIMER_RIGHT_0_TM'] * 10) / 10}"
    return row


def run_primer3(mut_df, chroms_folder='', pcr_config=PCR_config, primer3_config=primer3_config):
    '''
    input is df with columns 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene' with optional id columns to the left
    output is id columns + 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene' + primer cols
    primer cols:
        
    '''
    
    # apply pcr size to primer3_config
    primer3_config['PRIMER_PRODUCT_SIZE_RANGE'] = [pcr_config['prod_size_min'], pcr_config['prod_size_max']]
    primer3_config.update(pcr_config)

    # adjust column dtypes
    mut_df.loc[:, 'Chr'] = mut_df['Chr'].astype('str')
    mut_df = mut_df.loc[:,:'Alt']
    org_cols = list(mut_df.columns)
    df_list = []
    # cycle through (formatted) chromosomes
    # + load chromosome sequence
    # + create primer_df for mutations on that chromosome
    # + concat all mutations 
    for chrom in mut_df['Chr'].unique():
        chromo = get_chrom(chrom, chroms_folder=chroms_folder)
        # select the mutations on current chrom
        chr_df = mut_df.query('Chr == @chrom')
        print(f"Processing {len(chr_df.index)} mutations on {chrom}." )
        primer_df = chr_df.apply(get_primer_pair, chrom_dict=chromo, config=primer3_config, axis=1)
        df_list.append(primer_df)
    primer_df = pd.concat(df_list, sort=True)
    primer_df = primer_df[org_cols + ['fwdPrimer', 'revPrimer', 'Status', 'Temp', 'AmpliconRange', 'AmpliconSize', 'InsertRange', 'InsertSize', 'InsertSeq', 'offsetL', 'offsetR']]
    return primer_df
