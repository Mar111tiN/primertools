import pandas as pd
from multiprocessing import Pool
from functools import partial

from p3_core import compute_primers
from genomics import get_chrom
from script_utils import show_output

# ############ DEFAULT CONFIGS ##################
# for description of config see:
#       manual in https://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_manual.htm#PRIMER_MAX_SELF_ANY
PCR_config = {
    "seq_len": 500,
    "center_offSet": 5,
    "prod_size_min": 120,
    "prod_size_max": 220
    }

p3_default_config = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_INTERNAL_MAX_SELF_END": 8,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 55.0,
    "PRIMER_MAX_TM": 65.0,
    "PRIMER_MIN_GC": 20.0,
    "PRIMER_MAX_GC": 80.0,
    "PRIMER_MAX_POLY_X": 100,
    "PRIMER_INTERNAL_MAX_POLY_X": 100,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_MAX_SELF_ANY": 12,
    "PRIMER_MAX_SELF_END": 8,
    "PRIMER_PAIR_MAX_COMPL_ANY": 12,
    "PRIMER_PAIR_MAX_COMPL_END": 8,
}

p3_improved_config = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 55.0,
    "PRIMER_MAX_TM": 65.0,
    "PRIMER_MIN_GC": 20.0,
    "PRIMER_MAX_GC": 80.0,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_INTERNAL_MAX_SELF_END": 8,
    "PRIMER_MAX_POLY_X": 100,
    "PRIMER_INTERNAL_MAX_POLY_X": 100,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    # set 1 to actually use the thermodynamic calculations
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    "PRIMER_MAX_SELF_ANY": 12,
    "PRIMER_MAX_SELF_END": 8,
    "PRIMER_WT_SELF_END": 1,  # use Primer_max_self_end
    "PRIMER_MAX_SELF_END_TH": 30,
    "PRIMER_WT_SELF_END_TH": 1,  # Primer_max_self_end_th
    "PRIMER_PAIR_MAX_COMPL_ANY": 12,
    "PRIMER_PAIR_MAX_COMPL_END": 8,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8,
    "PRIMER_PAIR_WT_COMPL_ANY": 2,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 30,
    "PRIMER_PAIR_WT_COMPL_ANY_TH": 2,
    "PRIMER_MAX_HAIRPIN_TH": 47,
    "PRIMER_WT_HAIRPIN_TH": 1,
}


def get_primer_df(mut_df, primer3_config, chroms_folder, chrom):
    """
    allocates dfs chrom-wise and controls row-wise computation
    """
    show_output(f"Running primer3 for chromosome {chrom}.", multi=True)
    chrom_dict = get_chrom(chrom, chroms_folder)
    chr_df = mut_df.query("Chr == @chrom")
    primer_df = chr_df.apply(
        compute_primers, chrom=chrom_dict, config=primer3_config, axis=1
    )
    show_output(f"Finished chrom {chrom}.", multi=True)
    return primer_df


def run_primer3(
    mut_df,
    chroms_folder=".",
    pcr_config={},  # use defaults defined at top
    primer3_config={},  # use defaults defined at top
    threads=1,
):
    """
    input is df with columns 'Chr', 'Start', 'End', 'Ref', 'Alt' with optional id columns (everything left of Chr)
    output is id columns + 'Chr', 'Start', 'End', 'Ref', 'Alt' + primer cols
    primer cols:

    """
    
    # apply pcr size to primer3_config
    primer3_config["PRIMER_PRODUCT_SIZE_RANGE"] = [
        pcr_config["prod_size_min"],
        pcr_config["prod_size_max"],
    ]
    primer3_config.update(pcr_config)

    mut_df.loc[:, "Chr"] = mut_df["Chr"].astype("str")

    # COLS
    base_cols = ["Chr", "Start", "End", "Ref", "Alt"]
    # keep possible columns left of Chr
    # save the id columns into org_df for later merge
    keep_cols = list(mut_df.columns[: list(mut_df.columns).index("Chr")])

    keep_df = mut_df.loc[:, keep_cols + base_cols]

    mut_df = mut_df.loc[:, base_cols]

    df_list = []
    # cycle through (formatted) chromosomes
    # + load chromosome sequence
    # + create primer_df for mutations on that chromosome
    # + concat all mutations

    # ##### MULTIPROCESSING
    chrom_list = mut_df["Chr"].unique()
    show_output(f"Allocating processor pool for {threads} threads.")
    pool = Pool(threads)
    df_list = pool.map(
        partial(get_primer_df, mut_df, primer3_config, chroms_folder), chrom_list
    )

    primer_df = pd.concat(df_list, sort=True)
    new_cols = [
        "fwdPrimer",
        "revPrimer",
        "Status",
        "Temp",
        "AmpliconRange",
        "AmpliconSize",
        "InsertRange",
        "InsertSize",
        "InsertSeq",
        "offsetL",
        "offsetR",
    ]
    for col in ["AmpliconSize", "InsertSize", "offsetL", "offsetR"]:
        primer_df[col] = primer_df[col].fillna(0).astype(int)

    primer_df = keep_df.merge(primer_df[base_cols + new_cols])

    return primer_df


def primer3_master(
    i,
    o,
    chroms_folder,
    threads,
    PCR_config=PCR_config,
    primer3_config=p3_improved_config,
):
    """
    wrapper around run_primer3 that allows for injecting with adjusted PCR and Primer3_configs
    and controls input and output
    """

    # #### LOAD file ###################
    show_output(f"Loading {i} for primer3 computation. ", end="")
    filter1_df = pd.read_csv(i, sep="\t")
    show_output(f"{len(filter1_df.index)} mutations found.", time=False)

    # #### run primer3 #################
    primer_df = run_primer3(
        filter1_df,
        chroms_folder,
        pcr_config=PCR_config,
        primer3_config=primer3_config,
        threads=threads,
    )

    # ###### write to file #############
    primer_df.drop_duplicates().to_csv(o, sep="\t", index=False)
    show_output(f"Writing primer list to {o}.")
