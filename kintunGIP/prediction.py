from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import repeatfinder as rf
import pandas as pd
from ast import literal_eval
import numpy as np
import collections
import pyhmmer
from tqdm import tqdm


def extract_cds_from_range(genbank_file, start, end):
    cds_records = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Check if the CDS feature overlaps with the specified range
                if start <= feature.location.end and end >= feature.location.start:
                    # Extract the CDS sequence and translate it to protein
                    cds_seq = feature.extract(record.seq)
                    protein_seq = cds_seq.translate(table=11,
                                                    to_stop=True)  # Translate to protein, table=11 is for Bacterial
                    if len(protein_seq) > 0:
                        cds_records.append(
                            SeqRecord(protein_seq, id=feature.qualifiers.get("protein_id", ["unknown_protein"])[0]))

    return cds_records


def write_protein_fasta(cds_records):
    with open("tmp_file.fasta", "w") as temp_fasta:
        SeqIO.write(cds_records, temp_fasta, "fasta")
        temp_fasta_name = temp_fasta.name
    return temp_fasta_name


# here function that returns hits
def check_integrases_presence(genbank_file, start, end, hmm):
    cds_records = extract_cds_from_range(genbank_file, start, end)
    if len(cds_records) == 0:
        return []
    else:
        temp_fasta_file = write_protein_fasta(cds_records)
        with pyhmmer.easel.SequenceFile(temp_fasta_file, digital=True) as seq_file:
            sequences = seq_file.read_block()
        pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
        hits = pipeline.search_hmm(hmm, sequences)
        Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
        results = []
        for hits in pyhmmer.hmmsearch(hmm, sequences):
            cog = hits.query_name.decode()
            for hit in hits:
                if hit.included:
                    results.append(Result(hit.name.decode(), cog, hit.score))
        best_results = {}
        keep_query = set()
        for result in results:
            if result.query in best_results:
                previous_bitscore = best_results[result.query].bitscore
                if result.bitscore > previous_bitscore:
                    best_results[result.query] = result
                    keep_query.add(result.query)
                elif result.bitscore == previous_bitscore:
                    if best_results[result.query].cog != hit.cog:
                        keep_query.remove(result.query)
            else:
                best_results[result.query] = result
                keep_query.add(result.query)
        filtered_results = [best_results[k] for k in sorted(best_results) if k in keep_query]
        if len(filtered_results) > 0:
            return filtered_results
        else:
            return []


def check_overlap(range1, range2):
    """Check if two ranges overlap."""
    # If range1 starts after range2 ends or range1 ends before range2 starts, they do not overlap
    if range1[1] < range2[0] or range2[1] < range1[0]:
        return None  # No overlap, return None
    else:
        return (max(range1[0], range2[0]), min(range1[1], range2[1]))  # Return the overlapping range


def check_tuple_overlap(tuple_values, range_to_check):
    """Check if the ranges delimited by the first and second or the third and fourth values of a tuple overlap with a given range."""
    # Check if the range delimited by the first and second values overlaps with the given range
    overlap1 = check_overlap((tuple_values[0], tuple_values[1]), range_to_check)
    if overlap1:
        return tuple_values  # Return the original tuple if there's an overlap with the first range
    # Check if the range delimited by the third and fourth values overlaps with the given range
    overlap2 = check_overlap((tuple_values[2], tuple_values[3]), range_to_check)
    if overlap2:
        return tuple_values  # Return the original tuple if there's an overlap with the second range
    return None  # No overlap found, return None


def check_repeats_presence(genbank_file, start, end, sense, locus, gene_size):
    record = [record for record in SeqIO.parse(open(genbank_file), "gb")][0]
    # ----(tDNA)---------------->
    if (sense == 1 and locus == "DR") or (sense == -1 and locus == "UR"):
        new_start = start - 50 - gene_size
        new_end = end + 50
        range_start = 0
        range_end = gene_size + 100

    # <---------------------(tDNA)----
    elif (sense == 1 and locus == "UR") or (sense == -1 and locus == "DR"):
        new_start = start - 50
        new_end = end + 50 + gene_size
        range_start = (new_end - new_start) - 100 - gene_size
        range_end = (new_end - new_start)

    repeats = rf.get_repeats(str(record.seq[new_start:new_end]), gap=3)
    filtered_repeats = [check_tuple_overlap(tuple_values, (range_start, range_end)) for tuple_values in repeats]

    filtered_repeats = [i for i in filtered_repeats if i is not None]
    filtered_repeats = [tuple_values for tuple_values in filtered_repeats if
                        tuple_values[1] > tuple_values[0] and tuple_values[3] > tuple_values[2]]
    filtered_repeats = [tuple_values for tuple_values in filtered_repeats if
                        max([tuple_values[2] - tuple_values[1], tuple_values[3] - tuple_values[0]]) > 3000]
    filtered_repeats = [tuple_values for tuple_values in filtered_repeats if
                        max([tuple_values[1] - tuple_values[0], tuple_values[3] - tuple_values[2]]) > 15]

    if len(filtered_repeats) > 0:
        to_return = []
        for tuple_values in filtered_repeats:
            list_return = []
            for coord in tuple_values:
                list_return.append(coord + new_start)
            to_return.append(tuple(list_return))
        return to_return
    else:
        return []


def check_gc_content(genbank_file, start, end):
    record = [record for record in SeqIO.parse(open(genbank_file), "gb")][0]
    pilr = record[start:end]
    if gc_fraction(pilr.seq) < gc_fraction(record.seq) * 0.95:
        return "Lower"
    else:
        return "Not Lower"


def check_properties(row, sense, indir, ext, hmm):
    if row.Prevalence > 0.5:
        if sense == "DR" and not pd.isna(row.PILR_DR):
            start, end = literal_eval(row.PILR_DR)
        elif sense == "UR" and not pd.isna(row.PILR_UR):
            start, end = literal_eval(row.PILR_UR)
        else:
            return ()
        genbank_file = f"{indir}/{row.File}.{ext}"
        check_int = check_integrases_presence(genbank_file, start, end, hmm)
        check_rep = check_repeats_presence(genbank_file, start, end, row.Strand, sense, abs(row.End - row.Start))
        check_gcc = check_gc_content(genbank_file, start, end)
        if check_int == [] and check_rep == []:
            return ()
        else:
            return (check_int, check_rep, check_gcc)
    else:
        return ()


def prediction(indir, output_folder, prefix, extension):
    tqdm.pandas()
    tdnas_df = pd.read_csv(f"{output_folder}/{prefix}_tDNAs.csv", header=0, index_col=0)
    cds_df = pd.read_csv(f"{output_folder}/{prefix}_core_CDSs.csv", header=0, index_col=0)
    panaroo_roary = pd.read_csv(f"{output_folder}/{prefix}_coregenome_detail.csv", header=0, index_col=0)
    with pyhmmer.plan7.HMMFile("int_bact_pha.hmm") as hmm_file:
        hmm = hmm_file.read()
    tdnas_df["diag_UR"] = tdnas_df.progress_apply(check_properties, args=("UR", indir, extension, hmm,), axis=1)
    tdnas_df["diag_DR"] = tdnas_df.progress_apply(check_properties, args=("DR", indir, extension, hmm,), axis=1)
    tdnas_df.to_csv(f"{output_folder}/{prefix}_EMGs_prediction.csv")

