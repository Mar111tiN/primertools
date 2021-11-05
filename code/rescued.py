#########################################
########  STUFF FROM OLD NOTEBOOKS  #####


def num2well(num):
    '''
    loads PCR setup for association sample :: MutID(PrimerName) -> SampleMutID
    '''
    
    string = "ABCDEFGH"
    char = string[(num - 1) % 8]
    count = int((num - 1) / 8) + 1
    return f"{char}{count}"


##### INDEXING #####
# ### check overlaps
# + illegal overlaps are between different samples at overlapping positions in the same index
# + df has to be sorted by index, amplicon start site and then by SampleID


def showOverlap(df, padding=10):
    '''
    checks if the df has any overlaps as described above
    '''

    full_df = df.copy()
    # resort
    full_df = full_df.sort_values(['index', "Chr", "AmpStart", "SampleID"])
    # make an integer for any kind of overlap
    full_df['internal_overlap'] = ((full_df['AmpStart'] + padding < full_df.shift(1)['AmpEnd']) & (full_df['Chr'] == full_df.shift(1)['Chr']) & (full_df['index'] == full_df.shift(1)['index']) & (full_df['SampleID'] != full_df.shift(1)['SampleID'])).astype(int)
    return full_df.reset_index(drop=True)

def hasOverlap(df, padding=10):
    '''
    boolean checks if the df has any overlaps as described above
    '''
    
    full_df = showOverlap(df)
    return len(full_df.query("internal_overlap == 1")) > 0


black_list = [f"NHL{i}-{t}" for i in [1,2,8] for t in ["WB", "GR"]]


def get_index_df(org_df):
    '''
    returns the sum of sequenczing samples per index
    '''
    
    # assign microliters
    df = org_df.copy()
    df['µl'] = 2
    df.loc[df['status'].str.startswith("F"), 'µl'] = 0
    df.loc[df['status'].str.startswith("P"), 'µl'] = 4
    index_df = df.loc[~df['status'].str.startswith("F")].groupby("index").agg({"MutID":"count", "Project":"first", "µl":"sum", "SampleID":"first"}).rename(dict(MutID="sampleCount"), axis=1).reset_index()
    
    return index_df


def get_min_index(df, black_list=black_list, offset=0):
    '''
    get the index containing fewest samples from samples exluding black_list
    '''
    # remove -1 index
    index_df = get_index_df(df.query("index != -1")).query("SampleID not in @black_list").sort_values("sampleCount", ascending=True).iloc[offset]['index']
    
    return index_df


### FROM VALIDATION

def get_internal_overlap(df, padding=30):  
    df = df.sort_values(['Chr', 'Start', 'Plate', 'Well']).reset_index(drop=True)
    df[['WellAH', 'Well19']] = df['Well'].str.extract(r'([A-H])([0-9]+)')
    df['Well19'] = df['Well19'].astype(int)
    
    # extract the positions from InsertRange
    df[['InsertChr', 'InsertStart', 'InsertEnd']] = df['InsertRange'].str.extract('(?P<Chr>chr[0-9XY]+):(?P<Start>[0-9]+)-(?P<End>[0-9]+)')
    for col in ['InsertStart', 'InsertEnd']:
        df[col] = df[col].astype(int)
        
    df['internal_overlap'] = ((df['InsertStart'] + padding < df.shift(1)['InsertEnd']) & (df['Chr'] == df.shift(1)['Chr'])).astype(int)
    df.loc[df['internal_overlap'] == 1,'ov_with'] = df.shift(1)['MutID']
    df['ov_with'] = df['ov_with'].fillna('')
    # expand the overlap to the row below
    df.loc[df['internal_overlap'] == 0, 'internal_overlap'] = df.shift(-1)['internal_overlap'].fillna(0).astype(int)
    # df['ov'] = df[((df['InsertStart'] + padding < df.shift(1)['InsertEnd']) & (df['Chr'] == df.shift(1)['Chr']))]
    return df


def get_overlap_groups(df):
    indexCount = 1
    index_df = pd.DataFrame()
    # for i, row in df.iterrows():
        # indexCount = row['internal_overlap'] * (indexCount + row['internal_overlap'])
        # row['indexCount'] = int(indexCount)
    has_index = df.query('internal_overlap == 1')
    for i, row in has_index.iterrows():
        if row['ov_with'] == '':
            indexCount += 1
            index = 1
        else:
            index += 1
        row['indexCount'] = index
        row['overlapGroup'] = indexCount
        index_df = index_df.append(row)
        
        
    # merge index info into df
    df = df.merge(index_df, how='left')
    df['overlapGroup'] = df['overlapGroup'].fillna(0).astype(int)
    df['indexCount'] = df['indexCount'].fillna(0).astype(int)
    return df # .sort_values(['Plate', 'WellAH', 'Well19']).reset_index(drop=True)


def get_primer_reuse(df):
    # prepare the df
    df = df.sort_values(['fwd-Primer', 'rev_Primer', 'Plate', 'WellAH', 'Well19']).reset_index(drop=True)
    
    df['reuse'] = df[(df['fwd-Primer'] != df.shift(1)['fwd-Primer']) | (df['rev_Primer'] != df.shift(1)['rev_Primer'])]['MutID']
    df['reuse'] = df['reuse'].fillna(method='ffill')
    df.loc[df['reuse'] == df['MutID'],'reuse'] = ''
    return df.sort_values(['Plate', 'WellAH', 'Well19'])


def get_index(df, padding=30):
    '''
    get the overlapping samples and the coords
    '''

    overlap_df = get_internal_overlap(df, padding=padding)
    group_df = get_overlap_groups(overlap_df)
    # increment the index for overlapping positions
    index_max = group_df['indexCount'].max()
    
    # spread the overlap mutations evenly across the indices
    has_index = group_df.query('indexCount > 0').reset_index(drop=True)
    has_index['LibIndex'] = has_index.index % index_max
    
    # divide the non-overlapping mutations in even chunks and apply LibIndex
    no_index = group_df.query('indexCount == 0').sort_values(['Plate', 'WellAH', 'Well19']).reset_index(drop=True)
    no_index['LibIndex'] = (no_index.index / (len(no_index.index) / index_max)).astype(int)
    df = pd.concat([has_index,no_index]).sort_values(['Plate', 'WellAH', 'Well19']).reset_index(drop=True)
    print(f"Found {len(df.query('internal_overlap == 1').index)} overlaps - requiring {index_max} indices")
    
    print('Making plate schemas')
    # make plate_dfs for the indices
    index_dict = get_plates(df, 'LibIndex')
    
    # analyze the the possible reuse of primers
    df = get_primer_reuse(df)
    primer_dict = get_plates(df, 'reuse')
    
    return df.drop(columns=['WellAH', 'Well19']), index_dict, primer_dict



# combined code 
# see NEBNextNeu.ipynb

def merge_indices(main_df, df, black_list=[]):
    '''
    takes an indexed main_df and tries to merge df into it
    1) non-overlapping rows in df are bundled into one index
    2) internally overlapping mutations are assigned to index in main_df with fewest samples
    '''
    # remove any duplicate SampleMutID
    # they do exist!! need to be seen by pipeline but not by indexing
    df = df.drop_duplicates("SampleMutID")

    # set arbitrary index for checking overlaps
    df['index'] = -1
    # assign internal overlaps
    df = showOverlap(df)
    
    #### ASSIGN NON-OVERLAPPING MUTATIONS TO NEW INDEX ########
    
    
    # make new index for non_overlapping samples
    next_index = get_index_df(main_df)['index'].max() + 1
    print("Assigning non-overlapping bulk to index = ", next_index)
    df.loc[df["internal_overlap"] == 0, "index"] = next_index
    
    #### ASSIGN OVERLAPPING MUTATIONS ########
    # add overlap_df to main_df and 
    all_df = pd.concat([main_df, df.query('internal_overlap == 1')]).reset_index(drop=True).drop('internal_overlap', axis=1)

    overlap_df = all_df.query('index == -1')
    print(f"Assigning overlapping {len(overlap_df.index)} samples to existing index groups")
    
    # step through rows
    for i, row in overlap_df.iterrows():
        # get the current index
        index = row.name
        new_index = get_min_index(all_df, black_list=black_list)
        print(f"Assigning SampleMutID {row['SampleMutID']} (row {index})")
        all_df.loc[index,'index'] = new_index
        offset = 1
        while hasOverlap(all_df.query("index > -1")):
            print(f"SampleMutID {row['SampleMutID']} hasOverlap in index {new_index}!")
            # move to next index in line using offset
            new_index = get_min_index(all_df, offset=offset, black_list=black_list)
            all_df.loc[index,'index'] = new_index
            offset += 1
        print(f"SampleMutID {row['SampleMutID']} has been assigned to index {new_index}!")
    
    # bringing in the non_overlapping bulk mutations
    total_df = pd.concat([all_df, df.query('index > -1').drop('internal_overlap', axis=1)]).reset_index(drop=True).sort_values(["index", "Well"])
    return total_df
    




def make_index_table(org_df):
    '''
    creates a table with a column per index
    '''
    
    df = add_info(org_df)
    # get the maximum rows for initiating a clean df
    max_rows = get_index_df(df)['sampleCount'].max()

    table_df = pd.DataFrame(index=range(max_rows))
    for i in df['index'].unique():
        table_df[f"index{i}"] = df.query('index == @i').reset_index()['info']
    return table_df