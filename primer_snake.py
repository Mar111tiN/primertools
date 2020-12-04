from code.p3 import primer3_master


def main(s):
    # ############ SNAKEMAKE ##################
    i = s.input
    o = s.output
    c = s.config
    cc = c["primer3"]

    primer3_master(
        s.input[0],
        s.output[0],
        chroms_folder=s.params.genome_split,
        threads=cc["threads"],
        PCR_config=cc["params"]
        # primer3_config taken from p3_improved_defaults
    )


if __name__ == "__main__":
    main(snakemake)
