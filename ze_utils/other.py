def identify_loops(pose):
    """
    Returns a list of loops, where each loop is a list of residue indexes.
    CAUTION: Doesn't work perfectly. Loop recognition is performed by pyrosetta
    in 'set_ss_from_phipsi', which isn't too good.
    """

    from pyrosetta.rosetta.core.fragment import SecondaryStructure
    from pyrosetta.rosetta.core.pose import set_ss_from_phipsi

    set_ss_from_phipsi(pose)

    full_ss = SecondaryStructure(pose)
    current_loop, raw_loops = [], []
    last_ss = 'X'
    for residue_index in range(1, len(pose.sequence())):
        ss = full_ss.secstruct(residue_index)
        if ss != last_ss:
            if last_ss == 'L':
                raw_loops.append(current_loop)
            current_loop = []
        if ss == 'L':
            current_loop.append(residue_index)
        last_ss = ss

    loops = []
    for loop in raw_loops:
        if len(loop) > 1:
            loops.append(loop)

    return loops