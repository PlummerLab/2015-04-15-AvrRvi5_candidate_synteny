import re

def read_promer_coords(coords_handle):
    """ Parse promer coords file.

    Keyword arguments:
    coords_file -- Path to promer output coords file (string, required)

    yields:
        label -- An integer, in ascending order of when they are encountered.
        psim -- % similarity of the alignment (based on the scoring matrix
            that you used in promer).
        pid -- % AA identity in the alignment.
        pstp -- % stop codons in the alignment
        <reference_id> -- A dictionary containing the seqid, start position,
            end position, and strand of the alignment for the reference
            sequence provided to promer.
        <query_id> -- As with 'reference' but for the promer query sequence.
        reference -- <reference_id>
        query -- <query_id>
    """

    start_finder = re.compile(r"=+")
    line_split = re.compile(r"\s+\|\s+|\s+")

    started = False
    for i, line in enumerate(coords_handle):
        if i == 0:
            genomes = line.split()
        line = line.strip()
        if not started:
            if start_finder.match(line) != None:
                started = True
            continue

        comp = dict()
        line = line_split.split(line)
        comp['label'] = i

        comp['pid'] = float(line[6])  # %identity
        comp['psim'] = float(line[7])  # %similarity
        comp['pstp'] = float(line[8])  # %stop codons

        comp[line[11]] = {
            "start": int(line[0]),
            "length": int(line[1]) - int(line[0]),
            "strand": int(line[9]),
            "seqid": line[11],
            }
        comp['reference'] = line[11]
        comp[line[12]] = {
            "start": int(line[2]),
            "length": int(line[3]) - int(line[2]),
            "strand": int(line[10]),
            "seqid": line[12]
            }
        comp['query'] = line[12]

        # Check to see if it is a self match
        if line[11] == line[12]:
            if (min(line[0], line[1]) == min(line[2], line[3])) and \
                    (max(line[0], line[1]) ==  max(line[2], line[3])):
                continue
        yield comp
    return
