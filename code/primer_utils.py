def mut2insert(mut={}, seq={}, return_none=False):
    '''
    takes a mutation dictionary of shape:
    {
        'Chrom': 'chr4,
        'Start': 13331,
        'End': 13331,
        'Ref': '-',
        'Alt': 'A'
        }
    and a seq dictionary of shape:
    {
        'Chrom': 'chr4,
        'Start': 13300,
        'End': 123123,
        'seq': 'ATTTCTCCCACTCCCACA',
        }
    and returns the seq with the mutation inserted into the sequence
    if mutation location is out of bounds of sequence, only the sequence is returned without editing
    '''

    if mut['Chr'] != seq['Chr']:
        return None if return_none else seq['seq']
    if (mut['Start'] < seq['Start']) or (mut['End'] > seq['End']):
        return None if return_none else seq['seq']
    bases = ['A', 'C', 'G', 'T']
    # case SNP
    start = mut['Start'] - seq['Start']
    end = mut['End'] - seq['Start']

    if mut['Ref'] == "-":
        upstream = seq['seq'][:start]
        downstream = seq['seq'][end:]
        mutation = f"<<+{mut['Alt']}>>"
    else:
        upstream = seq['seq'][:start]
        downstream = seq['seq'][end + 1:]
        if mut['Alt'] == "-":
            mutation = f">∆{mut['Ref']}<"
        else:
            # SNP
            mutation = f"({mut['Ref']}>{mut['Alt']})"
    return f"{upstream}{mutation}{downstream}"

    
def edit_seq(row):
    # convert mutation to dict
    def mut2dict(row):
        # mutation to dict
        chrom = row['Chr']
        start = row['Start']
        end = row['End']
        mut_dict = {
            'Chr': row['Chr'],
            'Start': row['Start'],
            'End': row['End'],
            'Ref': row['Ref'],
            'Alt': row['Alt']
        }
        return mut_dict
    mut_dict = mut2dict(row)
    
    insert_dict = {
        'Chr': row['Chr'],
        'Start': row['InsertStart'],
        'End': row['InsertEnd'],
        'seq': row['InsertSeq']
    }
    
    edited_seq = mut2insert(mut=mut_dict, seq=insert_dict, return_none=True)
    return edited_seq
