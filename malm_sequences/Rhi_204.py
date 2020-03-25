import numpy as np

def mk_seq(a_injection = 1, b_return = 64):
    """
    write the rhizotron malm sequence with 204 measurements
    parameters:
    ===
    * a_injection: number of the electrode used for injection (t: int)
    * b_return: number of the electrode used as return (t: int)
    returns:
    ===
    * seq: malm sequence (t: np.ndarray)
    """
    # check some args
    types_a_injection = [int]
    types_b_return = [int]
    if not type(a_injection) in types_a_injection:
        raise TypeError('inappropriate arg type, a_injection')
    if not type(b_return) in types_b_return:
        raise TypeError('inappropriate arg type, b_return')

    # function    
    a = a_injection
    b = b_return
    
    seq = []
    
    # horizontal
    electrodes = np.arange(1,65)
    electrodes = electrodes.reshape((8,8))
    for row in electrodes:
        for i in range(len(row) -1):
            m = row[i]
            n = row[i + 1]
            seq.append([a, b, m, n])

    # vertical
    rotated_electrodes = np.flipud(np.rot90(electrodes))
    for row in rotated_electrodes:
        for i in range(len(row) -1):
            m = row[i]
            n = row[i + 1]
            seq.append([a, b, m, n])

    # diag \ -6 : 6
    for d in range(-6, 7):
        diagonal = np.diag(electrodes, k = d)
        for i in range(len(diagonal) - 1):
            m = diagonal[i]
            n = diagonal[i + 1]
            seq.append([a, b, m, n])

    # diag / 
    flipped_electrodes = np.fliplr(electrodes)
    for d in range(-6, 7):
        diagonal = np.diag(flipped_electrodes, k = d)
        for i in range(len(diagonal) - 1):
            m = diagonal[i]
            n = diagonal[i + 1]
            seq.append([a, b, m, n])

    # remove where the same electrode appears twice
    sequence = np.array(seq)
    indeces_delete = []
    for i, row in enumerate(sequence):
        if len(np.unique(row)) != len(row):
            indeces_delete.append(i)
    
    sequence_clean = np.delete(sequence, indeces_delete, axis = 0)

    return(sequence_clean)

def mk_full_malm(seq, VRTe = range(65, 371)):
    """
    tile a sequence and sobstitute the injection electrode
    
    parameters:
    ===
    * seq: the sequence (t: np.array)
    * VRTe: list of VRTe to use as injection (t: range, int, list, np.ndarray)

    return:
    ===
    * VRTe_seq : complete sequence for VRTe simulation (t: np.ndarray)
    """
    # check some args
    types_seq = [np.ndarray]
    types_VRTe = [range, list, np.ndarray, int]
    if not type(VRTe) in types_VRTe:
        raise TypeError('inappropriate arg type, VRTe')
    if not type(seq) in types_seq:
        raise TypeError('inappropriate arg type, seq')
    if isinstance(VRTe, int):
        VRTe = [VRTe]
    if len(VRTe) == 0:
        raise ValueError('inappropriate arg value, VRTe length must be > 0')
    # function

    num_VRTe = len(VRTe)
    VRTe_seq = np.tile(seq, (num_VRTe, 1))
    
    for i, vrte in enumerate(VRTe):
        VRTe_seq[range(i * len(seq), (i + 1) * len(seq)), 0] = vrte

    return(VRTe_seq)


if __name__ == "__main__":
    seq = mk_seq()
    VRTe_seq = mk_full_malm(seq)
    np.savetxt('seq.txt', seq, fmt= '%i %i %i %i')
    np.savetxt('VRTe_seq.txt', VRTe_seq, fmt= '%i %i %i %i')