#!/usr/bin/env python3
"""
Pure utility functions for Scaffold_builder.

Contains all stateless, side-effect-free functions used by the scaffolding
pipeline: sequence manipulation, alignment, coordinate parsing, gap/overlap
logic, FASTA I/O, and statistics.
"""

import os
import logging

VERSION = "2.2"

DEFAULTS = {
    "-a": 95,
    "-p": "Scaffold",
    "-i": 80,
    "-g": 5000,
    "-b": 0,
    "-r": "",
    "-q": "",
    "-t": 300,
}

IUPAC = {
    "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
    "AN": "N", "GN": "N", "CN": "N", "TN": "N", "NT": "N",
    "TY": "C", "CY": "C",
}


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence.

    Args:
        sequence (str): DNA sequence (A/T/C/G/N, any case).

    Returns:
        str: Reverse complement in uppercase.
    """
    table = str.maketrans("ATCGN", "TAGCN")
    return sequence.upper().translate(table)[::-1]


def fasta2hash(fasta_path):
    """Parse a FASTA file into a dict keyed by sequence ID.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        dict[str, str]: Sequence ID → uppercase sequence.
    """
    sequences = {}
    with open(fasta_path) as fh:
        for record in fh.read().split(">")[1:]:
            lines = record.split("\n")
            seq_id = lines[0].split()[0].replace("|", "")
            sequences[seq_id] = "".join(lines[1:]).upper()
    return sequences


def needleman_wunsch(seq_a, seq_b, min_identity):
    """Align two sequences with the Needleman-Wunsch algorithm.

    Args:
        seq_a (str): First sequence (current contig terminus).
        seq_b (str): Second sequence (next contig terminus).
        min_identity (float): Minimum percent identity to build consensus.

    Returns:
        list: [aligned_a, aligned_b, percent_identity, consensus_sequence]
    """
    def ambiguity_code(s1, s2):
        """Build IUPAC consensus from two gapped aligned sequences.

        Args:
            s1 (str): Aligned sequence 1 (may contain '-').
            s2 (str): Aligned sequence 2 (may contain '-').

        Returns:
            str: Consensus using IUPAC ambiguity codes.
        """
        consensus = ""
        length = len(s1)
        for i, (a, b) in enumerate(zip(s1, s2)):
            pair = [a, b]
            if i < length * 0.25:
                if a == "-":
                    continue
                elif "-" not in pair:
                    letters = sorted(set(pair))
                    consensus += IUPAC.get("".join(letters), "N") if len(letters) > 1 else letters[0]
                else:
                    consensus += a
            elif i > length * 0.75:
                if b == "-":
                    continue
                elif "-" not in pair:
                    letters = sorted(set(pair))
                    consensus += IUPAC.get("".join(letters), "N") if len(letters) > 1 else letters[0]
                else:
                    consensus += b
            else:
                if "-" in pair:
                    pair = sorted(set(c for c in pair if c != "-"))
                    if len(pair) > 1:
                        consensus += IUPAC.get("".join(pair), "N")
                    elif pair:
                        consensus += pair[0]
                else:
                    consensus += a
        return consensus

    def score_matrix():
        """Build and fill the DP scoring matrix.

        Returns:
            list[list[int]]: Filled (len_b+1) × (len_a+1) DP matrix.
        """
        matrix = [[0] * (len(seq_a) + 1) for _ in range(len(seq_b) + 1)]
        for col in range(1, len(seq_a) + 1):
            for row in range(1, len(seq_b) + 1):
                match = 1 if seq_b[row - 1] == seq_a[col - 1] else 0
                matrix[row][col] = max(
                    matrix[row - 1][col - 1] + match,
                    matrix[row][col - 1],
                    matrix[row - 1][col],
                )
        return matrix

    def traceback(matrix):
        """Trace back through the DP matrix to recover the alignment.

        Args:
            matrix (list[list[int]]): Filled DP scoring matrix.

        Returns:
            list: [aligned_a, aligned_b, percent_identity, consensus]
        """
        c, l = len(seq_a), len(seq_b)
        al_a = al_b = ""
        while l != 0 and c != 0:
            diag, front, up = matrix[l - 1][c - 1], matrix[l][c - 1], matrix[l - 1][c]
            best = max(diag, front, up)
            if best == diag:
                c -= 1; l -= 1
                al_a = seq_a[c] + al_a; al_b = seq_b[l] + al_b
            elif best == front:
                c -= 1
                al_a = seq_a[c] + al_a; al_b = "-" + al_b
            else:
                l -= 1
                al_a = "-" + al_a; al_b = seq_b[l] + al_b
            if l == 0 and c != 0:
                missing = len(seq_a) - len(al_a.replace("-", ""))
                al_a = seq_a[:missing] + al_a
                al_b = "-" * missing + al_b
                break
        matches = sum(a == b for a, b in zip(al_a, al_b))
        perc = matches * 100.0 / len(al_a) if al_a else 0.0
        consensus = ambiguity_code(al_a, al_b) if perc >= min_identity else ""
        return [al_a, al_b, perc, consensus]

    return traceback(score_matrix())


def coord2hash(coords_path, parameters, hash_sequences):
    """Parse a MUMmer .coords file into a coordinate hash.

    Sets parameters['begin_end'] if a contig maps at r1=1, splits it into
    _begin/_end entries, and stores a raw copy in parameters['rawMapping'].

    Args:
        coords_path (str): Path to the .coords file.
        parameters (dict): Updated in place: 'begin_end', 'rawMapping'.
        hash_sequences (dict): Sequence dict; updated if begin_end detected.

    Returns:
        dict[str, list]: Contig ID → list of [r1, r2, q1, q2] hit lists.
    """
    coord_hash = {}
    with open(coords_path) as fh:
        for line in fh:
            parts = line.replace("|", "").split()
            if len(parts) == 9:
                coords = [int(x) for x in parts[:4]]
                seq_id = parts[-1]
                if coords[0] == 1:
                    parameters["begin_end"] = seq_id
                coord_hash.setdefault(seq_id, []).append(coords)
    if parameters["begin_end"] != 0:
        coord_hash = _begin_end(coord_hash, parameters, hash_sequences)
    parameters["rawMapping"] = dict(coord_hash)
    return coord_hash


def _begin_end(coord_hash, parameters, hash_sequences):
    """Split a circular-reference contig into _begin and _end entries.

    Args:
        coord_hash (dict): Current coordinate hash.
        parameters (dict): Contains 'begin_end' key.
        hash_sequences (dict): Updated in place with _begin/_end sequences.

    Returns:
        dict: Updated coordinate hash.
    """
    be = parameters["begin_end"]
    mapping = coord_hash[be]

    def _extract(hit):
        _, _, q1, q2 = hit
        seq = hash_sequences[be]
        return reverse_complement(seq)[q2 - 1:q1] if q1 > q2 else seq[q1 - 1:q2]

    coord_hash[be + "_begin"] = [mapping[0]]
    hash_sequences[be + "_begin"] = _extract(mapping[0])
    coord_hash[be + "_end"] = [mapping[-1]]
    hash_sequences[be + "_end"] = _extract(mapping[-1])
    del coord_hash[be]
    return coord_hash


def clean_coords(coord_hash, parameters, hash_sequences):
    """Remove ambiguously mapped contigs and flatten single-hit entries.

    A contig is ambiguous if its second-best hit covers >= -a% of its length.

    Args:
        coord_hash (dict): Coordinate hash from coord2hash.
        parameters (dict): 'ambiguous' list updated in place.
        hash_sequences (dict): Used for contig length lookups.

    Returns:
        dict: Cleaned coordinate hash (each value is a flat [r1,r2,q1,q2]).
    """
    ambiguous = []
    for seq_id, hits in coord_hash.items():
        if len(hits) > 1 and seq_id != parameters["begin_end"]:
            sorted_hits = sorted(hits, key=lambda h: abs(h[2] - h[3]) + 1, reverse=True)
            second_len = abs(sorted_hits[1][2] - sorted_hits[1][3]) + 1
            if (second_len / len(hash_sequences[seq_id])) * 100 >= parameters["-a"]:
                parameters["ambiguous"].append(seq_id)
                ambiguous.append(seq_id)
            else:
                coord_hash[seq_id] = sorted_hits[0]
        else:
            coord_hash[seq_id] = hits[0]
    for seq_id in ambiguous:
        del coord_hash[seq_id]
    return coord_hash


def extend_mapping(coord_hash, hash_sequences):
    """Extend reference coordinates to account for unmapped contig termini.

    Args:
        coord_hash (dict): Cleaned coordinate hash (flat hits).
        hash_sequences (dict): Used for contig length lookups.

    Returns:
        dict: Updated coordinate hash with extended reference positions.
    """
    for seq_id, hit in coord_hash.items():
        if "_begin" in seq_id or "_end" in seq_id:
            continue
        contig_len = len(hash_sequences[seq_id])
        r1, r2, q1, q2 = hit
        if (q1 + q2 - 1) == contig_len:
            continue
        if q1 < q2:
            hit[0] = max(1, r1 - (q1 - 1))
            hit[1] = r2 + (contig_len - q2)
        else:
            hit[1] = r2 + (q2 - 1)
            hit[0] = max(1, r1 - (contig_len - q1))
        coord_hash[seq_id] = hit
    return coord_hash


def reverse_mapped(mapping, hash_sequences):
    """Reverse-complement contigs that mapped in reverse orientation (q1 > q2).

    Args:
        mapping (list): Sorted [[contig_id, hit], ...] pairs.
        hash_sequences (dict): Updated in place.
    """
    for seq_id, hit in mapping:
        if "_begin" not in seq_id and "_end" not in seq_id:
            if hit[2] > hit[3]:
                hash_sequences[seq_id] = reverse_complement(hash_sequences[seq_id])


def _remove_subset(sorted_list, parameters):
    """Remove contigs whose mapping lies entirely within another contig's region.

    Args:
        sorted_list (list): Sorted [[contig_id, hit], ...] pairs.
        parameters (dict): 'subSet' list updated in place.

    Returns:
        list: Filtered list.
    """
    c = 0
    while c < len(sorted_list) - 1:
        curr_hit, next_hit = sorted_list[c][1], sorted_list[c + 1][1]
        if curr_hit[0] < next_hit[0] and curr_hit[1] > next_hit[1]:
            parameters["subSet"].append([sorted_list[c + 1], sorted_list[c]])
            sorted_list.pop(c + 1)
            c -= 1
        c += 1
    return sorted_list


def write_fasta(seq_id, sequence, parameters):
    """Append a FASTA record to the scaffold output file.

    Args:
        seq_id (str): Sequence identifier.
        sequence (str): Nucleotide sequence.
        parameters (dict): Uses '-p' for output prefix.
    """
    with open(parameters["-p"] + "_Scaffold.fasta", "a") as fh:
        fh.write(">" + seq_id + "\n" + sequence.replace("\n", "") + "\n")


def write_overlap_alignment(id1, seq1, id2, seq2, parameters, action="Overlap_"):
    """Write a pairwise overlap or alignment FASTA file.

    Args:
        id1 (str): First sequence identifier.
        seq1 (str): First sequence (raw or aligned).
        id2 (str): Second sequence identifier.
        seq2 (str): Second sequence (raw or aligned).
        parameters (dict): Uses '-p' for directory prefix.
        action (str): File name prefix ('Overlap_' or 'Alignment_').
    """
    path = f"{parameters['-p']}_overlap_alignment/{action}{id1}_{id2}.fasta"
    with open(path, "w") as fh:
        fh.write(f">{id1}\n{seq1}\n>{id2}\n{seq2}\n")


def not_mapped(parameters, hash_sequences):
    """Identify unmapped contigs and write them to the scaffold FASTA.

    Args:
        parameters (dict): Contains 'rawMapping' and 'begin_end'.
        hash_sequences (dict): Full sequence dictionary.

    Returns:
        list[str]: IDs of contigs not present in rawMapping.
    """
    unmapped = []
    for seq_id in hash_sequences:
        if seq_id not in parameters["rawMapping"] and seq_id != parameters["begin_end"]:
            unmapped.append(seq_id)
            write_fasta(seq_id + "_not_mapped", hash_sequences[seq_id], parameters)
    return unmapped


def _fasta_stats(fasta_path):
    """Compute basic statistics for a FASTA file.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        tuple: (total_length, count, longest, all_lengths)
    """
    lengths = []
    with open(fasta_path) as fh:
        for record in fh.read().split(">")[1:]:
            lengths.append(len("".join(record.split("\n")[1:])))
    total = sum(lengths)
    return total, len(lengths), max(lengths) if lengths else 0, lengths


def _n50(lengths):
    """Compute the N50 value from a list of sequence lengths.

    Args:
        lengths (list[int]): Sequence lengths.

    Returns:
        int: N50 value.
    """
    expanded = sorted(length for length in lengths for _ in range(length))
    mid = len(expanded) // 2
    return (expanded[mid - 1] + expanded[mid]) // 2 if len(expanded) % 2 == 0 else expanded[mid]


def statistics(parameters):
    """Compute assembly vs scaffold statistics.

    Args:
        parameters (dict): Contains '-q' and '-p'.

    Returns:
        list: [asm_len, scaf_len, asm_count, scaf_count,
               asm_longest, scaf_longest, asm_n50, scaf_n50]
    """
    asm = _fasta_stats(parameters["-q"])
    scaf = _fasta_stats(parameters["-p"] + "_Scaffold.fasta")
    return [asm[0], scaf[0], asm[1], scaf[1], asm[2], scaf[2], _n50(asm[3]), _n50(scaf[3])]


def write_output(final_scaffold, parameters, hash_sequences):
    """Write the full scaffolding log and statistics to the output text file.

    Args:
        final_scaffold (list): Sorted [[contig_id, hit], ...] pairs.
        parameters (dict): Full parameters dict.
        hash_sequences (dict): Sequence dictionary.
    """
    info = f"Scaffold_builder version {VERSION} log file\n\nReference sequence: {parameters['-r']}\n\n"

    if parameters["begin_end"] != 0:
        info += (f"{parameters['begin_end']} maps to the beginning and end of the "
                 "reference, suggesting that this reference is circular.\n\n")

    if len(parameters["ambiguous"]) > 1:
        info += (f"Contigs mapped ambiguously: The contigs below map to the reference "
                 f"at least 2 times over >= {parameters['-a']}% of its length, "
                 "so it is not scaffolded.\n")
        for seq_id in parameters["ambiguous"]:
            write_fasta(seq_id + "_ambiguous_mapping", hash_sequences[seq_id], parameters)
            info += seq_id + "\n"

    if len(parameters["subSet"]) > 1:
        info += "\nContigs mapped entirely within a region where another contig has already been scaffolded are removed:\n\n"
        for pair in parameters["subSet"]:
            labels = [f"{p[0]} ({p[1][0]}-{p[1][1]})" for p in pair]
            subset_id = labels[0].split(" (")[0]
            write_fasta(subset_id + "_overlapping_hits_sub_region", hash_sequences[subset_id], parameters)
            info += f"{labels[0]} lies within the region of {labels[1]}\n"

    info += "\nFinal scaffolding:\n\n"
    info += "#contig\t\t\t\t\t5'ref\t\t\t\t\t3'ref\t\t\t\t\t5'contig\t\t\t\t\t3'contig\t\t\t\t\tlength\n"
    for entry in final_scaffold:
        hit = entry[1]
        info += "\t\t\t\t\t".join([entry[0]] + [str(x) for x in hit] + [str(abs(hit[-1] - hit[-2]) + 1)]) + "\n"

    info += "\nContigs that were not mapped to the reference by Nucmer\n"
    for seq_id in not_mapped(parameters, hash_sequences):
        info += seq_id + "\n"

    stats = statistics(parameters)
    info += "\nScaffolding statistics:\n\t\t\t\t\tAssembly\tScaffold\n"
    info += f"Total length\t{stats[0]}\t{stats[1]}\n"
    info += f"Number of sequences\t{stats[2]}\t{stats[3]}\n"
    info += f"Average length\t{round(stats[0]/stats[2], 2)}\t{round(stats[1]/stats[3], 2)}\n"
    info += f"Longest contig\t{stats[4]}\t{stats[5]}\n"
    info += f"N50\t{stats[6]}\t{stats[7]}\n"
    info += f"Non-overlapping contig pairs\t-\t{len(parameters['gaps'])}\n"
    info += f"Total length of gaps\t-\t{sum(parameters['gaps'])}\n"
    gaps = parameters["gaps"]
    if gaps:
        info += f"Longest gap\t-\t{max(gaps)}\n"
        info += f"Shortest gap\t-\t{min(gaps)}\n"
        info += f"Average gap size\t-\t{round(sum(gaps)/len(gaps), 2)}\n"
    else:
        info += "Longest gap\t-\t0\nShortest gap\t-\t0\nAverage gap size\t-\t0\n"
    info += f"Overlapping contig pairs (>=80% id)\t-\t{parameters['overlapOver']}\n"
    info += f"Overlapping contig pairs (<80% id)\t-\t{parameters['overlapBelow']}"

    with open(parameters["-p"] + "_output.txt", "w") as fh:
        fh.write(info)


def scaffold_sequences(mapping, parameters, hash_sequences):
    """Build scaffold sequences from the sorted, cleaned mapping.

    Fills gaps with Ns, resolves overlaps via Needleman-Wunsch, and writes
    scaffold FASTA records and overlap alignment files.

    Args:
        mapping (list): Sorted [[contig_id, hit], ...] pairs.
        parameters (dict): 'gaps', 'overlapOver', 'overlapBelow' updated in place.
        hash_sequences (dict): Updated in place for merged overlaps.
    """
    scaffold_seq = ""
    scaffold_num = 1
    next_element = None

    for c in range(len(mapping) - 1):
        curr, nxt = mapping[c], mapping[c + 1]
        next_element = nxt
        gap = nxt[1][0] - curr[1][1] - 1

        if gap >= 0:
            parameters["gaps"].append(gap)
            if gap >= parameters["-g"] and scaffold_seq:
                write_fasta("Scaffold_" + str(scaffold_num), scaffold_seq, parameters)
                scaffold_seq = ""
                scaffold_num += 1
            else:
                scaffold_seq += hash_sequences[curr[0]] + "N" * gap
        else:
            overlap_len = min(abs(gap), parameters["-t"])
            curr_seq, nxt_seq = hash_sequences[curr[0]], hash_sequences[nxt[0]]
            curr_overlap = curr_seq[-overlap_len:]
            nxt_overlap = nxt_seq[:overlap_len]
            curr_no_overlap = curr_seq[:-overlap_len]
            nxt_no_overlap = nxt_seq[overlap_len:]

            alignment = needleman_wunsch(curr_overlap, nxt_overlap, parameters["-i"])

            if alignment[2] < parameters["-i"] and overlap_len < parameters["-t"]:
                candidates = [
                    (needleman_wunsch(curr_seq[-ext:], nxt_seq[:ext], parameters["-i"]), ext)
                    for ext in range(overlap_len + 10, parameters["-t"] + 1, 10)
                ]
                if candidates:
                    best_al, best_ext = max(candidates, key=lambda x: x[0][2])
                    if best_al[2] > alignment[2]:
                        alignment = best_al
                        overlap_len = best_ext
                        curr_overlap = curr_seq[-overlap_len:]
                        nxt_overlap = nxt_seq[:overlap_len]
                        curr_no_overlap = curr_seq[:-overlap_len]
                        nxt_no_overlap = nxt_seq[overlap_len:]

            if alignment[2] >= parameters["-i"]:
                parameters["overlapOver"] += 1
                scaffold_seq += curr_no_overlap + alignment[3]
                hash_sequences[nxt[0]] = nxt_no_overlap
            else:
                parameters["overlapBelow"] += 1
                scaffold_seq += hash_sequences[curr[0]]
                if parameters["-b"] != 0 and scaffold_seq:
                    write_fasta("Scaffold_" + str(scaffold_num), scaffold_seq, parameters)
                    scaffold_seq = ""
                    scaffold_num += 1

            write_overlap_alignment(curr[0], curr_overlap, nxt[0], nxt_overlap, parameters)
            write_overlap_alignment(curr[0], alignment[0], nxt[0], alignment[1], parameters, "Alignment_")

    if next_element is not None:
        scaffold_seq += hash_sequences[next_element[0]]
    if scaffold_seq:
        write_fasta("Scaffold_" + str(scaffold_num), scaffold_seq, parameters)

    write_output(mapping, parameters, hash_sequences)


def sort_mapping(coord_hash, parameters, hash_sequences):
    """Sort, filter subsets, and scaffold the coordinate hash.

    Args:
        coord_hash (dict): Cleaned coordinate hash (flat hits).
        parameters (dict): Full parameters dict.
        hash_sequences (dict): Sequence dictionary.
    """
    coord_hash = extend_mapping(coord_hash, hash_sequences)
    sorted_list = sorted(
        ([seq_id, hit] for seq_id, hit in coord_hash.items()),
        key=lambda x: x[1][0],
    )
    reverse_mapped(sorted_list, hash_sequences)
    sorted_list = _remove_subset(sorted_list, parameters)
    scaffold_sequences(sorted_list, parameters, hash_sequences)
    logging.info("Scaffolding complete.")
