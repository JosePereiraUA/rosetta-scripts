# Questions: jose.manuel.pereira@ua.pt

from pymol import cmd

def get_sequence(query_object, verbose = True):
    """
    USAGE

    get_sequence OBJECT

    Returns a single string with the aminoacid sequence of the object.
    Similar to print(cmd.get_fastastr(OBJECT)) but does not separate and add
    chain ID names.
    """

    query_seq = ""

    # Retrieve the fasta string, while removing newlines
    fasta_seq = cmd.get_fastastr(query_object).replace("\n", "")

    # Remove the chain ID name headers
    for chain in cmd.get_chains(query_object):
        fasta_seq = fasta_seq.replace(">%s_%s" % (query_object, chain), "")

    # Print (by default) and return the sequence
    if verbose:
        print(fasta_seq)
    return fasta_seq


def color_by_sequence_conservation(query_object, reference_object):
    """
    USAGE

    color_by_sequence_conservation QUERY_OBJ, REFERENCE_OBJ

    Colors both the query and reference objects based on the aminoacid sequence
    conservation between the two objects. Residues are colored green if they
    share the same aminoacid in that position in both objects, otherwise are
    colored red. Prints and returns the percentage of sequence conserved.
    """

    # Get the sequences of the query and reference object
    query_seq     = get_sequence(query_object,     False)
    reference_seq = get_sequence(reference_object, False)

    # Sequences of query and reference objects must be of the same length.
    # No sequence alignment is performed.
    assert len(query_seq) == len(reference_seq), \
        "Objects do not have the same number of aminoacids."

    # Iterate over the two objects and retrieve the residue list. This makes
    # possible the comparison of two sequences where the residues are in the
    # same order but with different numbering.
    resi_list = {"query": [], "reference": []}
    cmd.iterate("%s and name CA" % (query_object), \
        "query.append(resi)", space = resi_list)
    cmd.iterate("%s and name CA" % (reference_object), \
        "reference.append(resi)", space = resi_list)

    # Iterate over the extracted sequences and compare each aminoacid
    conserved = 0
    for index in range(len(query_seq)):
        if query_seq[index] == reference_seq[index]:
            # In conserved aminoacid pairs (same aminoacid in both objects), use
            # green color green
            conserved += 1
            color = "green"
        else:
            # In non conserved aminoacids, use color red
            color = "red"

        # Color both objects with the defined color
        cmd.color(color, "resi %d and %s" % \
            (int(resi_list["query"][index]), query_object))
        cmd.color(color, "resi %d and %s" % \
            (int(resi_list["reference"][index]), reference_object))

    # Print and return the percentage of sequence conserved
    psc = (conserved/len(query_seq)) * 100
    print("Percentage of sequence conserved: %5.2f%%" % (psc))
    return psc


def color_by_sequence_blosum62_score(query_object, reference_object):
    """
    USAGE

    color_by_sequence_blosum62_score QUERY_OBJ, REFERENCE_OBJ

    Colors both the query and reference objects based on the aminoacid sequence
    BLOSUM62 score between the two objects. Residues are colored from a gradient
    between red and green where red residues have the lower BLOSUM62 and green
    residues the best. Prints and returns the BLOSUM62 score.
    """

    import sys
    from os.path import expanduser as expusr
    sys.path.insert(1, expusr("~/Desktop/scripts/"))
    from ze_utils.common import read_matrix_from_txt_file, blosum62

    # Get the sequences of the query and reference object
    query_seq     = get_sequence(query_object,     False)
    reference_seq = get_sequence(reference_object, False)

    # Sequences of query and reference objects must be of the same length.
    # No sequence alignment is performed.
    assert len(query_seq) == len(reference_seq), \
        "Objects do not have the same number of aminoacids."

    # Iterate over the two objects and retrieve the residue list. This makes
    # possible the comparison of two sequences where the residues are in the
    # same order but with different numbering.
    resi_list = {"query": [], "reference": []}
    cmd.iterate("%s and name CA" % (query_object), \
        "query.append(resi)", space = resi_list)
    cmd.iterate("%s and name CA" % (reference_object), \
        "reference.append(resi)", space = resi_list)

    # Read the BLOSUM62 matrix from the default location
    matrix, entries = read_matrix_from_txt_file()

    # Calculate the BLOSUM62 score between the given query and reference
    # sequences
    b62_scr, b62_res_scr = blosum62(matrix, entries, query_seq, reference_seq)
    
    # Create a Red-Green Gradient with 16 steps. This corresponds to each value
    # possible in the BLOSUM62 matrix (from -4 to 11). 
    colors = set_colors_from_RG_gradient(16)

    # Iterate over each residue score in the objects and color them on the
    # corresponding color from the generated Reg-Green Gradient
    for index, res_scr in enumerate(b62_res_scr):
        cmd.color("RG%d" % (res_scr + 4), "resi %d and %s" % \
            (int(resi_list["query"][index]), query_object))
        cmd.color("RG%d" % (res_scr + 4), "resi %d and %s" % \
            (int(resi_list["reference"][index]), reference_object))

    # Print and return the BLOSUM62 score
    print("BLOSUM62 score: %5d" % (b62_scr))
    return b62_scr


def set_colors_from_RG_gradient(n_steps):
    """
    USAGE

    set_colors_from_RG_gradient 16

    Create a Red-Green Gradient with 'n_steps'. Each color added to PyMOL will
    be named 'RGi', where i is the index, 0 being the red and n_steps the green
    (Ex: RG6).
    """

    assert type(n_steps) == int, "Number of steps should be an integer"

    step = 255 / n_steps
    for i in range(n_steps):
        r = 255 - step * i
        g = step * i
        cmd.set_color("RG%d" % (i), [r, g, 0])


# Export functions to PyMOL
cmd.extend("get_sequence", get_sequence)
cmd.extend("set_colors_from_RG_gradient", set_colors_from_RG_gradient)
cmd.extend("color_by_sequence_conservation", color_by_sequence_conservation)
cmd.extend("color_by_sequence_blosum62_score", color_by_sequence_blosum62_score)