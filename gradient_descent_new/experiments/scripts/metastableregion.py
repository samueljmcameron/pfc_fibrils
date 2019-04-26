import numpy as np

def find_metastable(R_fw,R_bw,E_fw,E_bw):

    # find all spots where two distinct phases exist (two different radii values)
    bool_1 = np.isclose(R_fw,R_bw,rtol=1e-3)
    find_j = np.where(bool_1==False)[0]

    if find_j.size > 0:   # if two distinct phases exist, then:

        # smallest gamma value that they both exist at
        jsmall = find_j[0]

        # largest gamma value that they both exist at
        jlarge = find_j[-1]

        # find the point where the fwd E becomes larger than the bkwd E
        j = (np.argmin(np.abs(E_fw[jsmall:jlarge+1]-E_bw[jsmall:jlarge+1]))
             +len(E_fw[:jsmall]))

    return jsmall,jlarge,j
