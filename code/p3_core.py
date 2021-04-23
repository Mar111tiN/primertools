import primer3
from primer_utils import mut2insert


def compute_primers(row, chrom, config):
    """
    row_based primer3 worker
    returns the best primer pair for a given position
    return value is [fwd_seq, fwd_tmp, rev_seq, rev_tmp, prod_size]
    active chromosome sequence is global variable chrom
    """

    # load sequence
    pos = row["Start"]
    half_seq = int(config["seq_len"] / 2)
    seq_start = pos - half_seq
    seq_end = pos + half_seq
    seq = chrom["sequence"][seq_start:seq_end]
    pad = int(config["center_offSet"] / 2)
    half_size = int(config["prod_size_min"] / 2)

    # calculate the target_range as offSet from center (half)
    offSet = half_size - 20 - pad
    target_start = half_seq - offSet
    target = [target_start, offSet * 2]
    setting = {
        "SEQUENCE_ID": "asdf",
        "SEQUENCE_TEMPLATE": seq,
        "SEQUENCE_TARGET": target,
    }
    primers = primer3.bindings.designPrimers(setting, config)

    # return '--' if nothing was found
    if primers["PRIMER_PAIR_NUM_RETURNED"] == 0:
        row["fwd_seq"] = row["fwd_tmp"] = row["rev_seq"] = row["rev_tmp"] = row[
            "prod_size"
        ] = "--"
        return row

    # # get chrom coords
    amp_start = seq_start + primers["PRIMER_LEFT_0"][0] + 1
    amp_end = seq_start + primers["PRIMER_RIGHT_0"][0] + 1

    row["AmpliconRange"] = f"{chrom['name']}:{amp_start}-{amp_end}"

    insert_start = amp_start + primers["PRIMER_LEFT_0"][1]
    insert_end = amp_end - primers["PRIMER_RIGHT_0"][1]
    row["InsertRange"] = f"{chrom['name']}:{insert_start}-{insert_end}"
    row["InsertSize"] = insert_end - insert_start
    insert_seq = chrom["sequence"][insert_start - 1 : insert_end]
    row["InsertSeq"] = insert_seq

    mut_dict = {
        "Chr": row["Chr"],
        "Start": row["Start"],
        "End": row["End"],
        "Ref": row["Ref"],
        "Alt": row["Alt"],
    }

    insert_dict = {
        "Chr": row["Chr"],
        "Start": insert_start,
        "End": insert_end,
        "seq": insert_seq,
    }

    row["InsertSeq"] = mut2insert(mut=mut_dict, seq=insert_dict, return_none=True)
    row["offsetL"] = row["Start"] - insert_start
    row["offsetR"] = insert_end - row["End"]

    row["fwdPrimer"] = primers["PRIMER_LEFT_0_SEQUENCE"]
    row["revPrimer"] = primers["PRIMER_RIGHT_0_SEQUENCE"]
    row["AmpliconSize"] = primers["PRIMER_RIGHT_0"][0] - primers["PRIMER_LEFT_0"][0] + 1
    row["Status"] = "not established"
    row[
        "Temp"
    ] = f"(fwd={int(primers['PRIMER_LEFT_0_TM'] * 10) / 10}|rev={int(primers['PRIMER_RIGHT_0_TM'] * 10) / 10})"
    return row
