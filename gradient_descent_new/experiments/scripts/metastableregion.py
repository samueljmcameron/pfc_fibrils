import numpy as np

def find_metastable(obsfwd,obsbkwd):

    # find all spots where two distinct phases exist (two different radii values)
    bool_1 = np.isclose(obsfwd.R(),obsbkwd.R(),rtol=1e-3)
    find_j = np.where(bool_1==False)[0]

    if find_j.size > 0:   # if two distinct phases exist, then:

        # smallest gamma value that they both exist at
        jsmall = find_j[0]

        # largest gamma value that they both exist at
        jlarge = find_j[-1]

        # find the point where the fwd E becomes larger than the bkwd E
        j = (np.argmin(np.abs(obsfwd.E()[jsmall:jlarge+1]-obsbkwd.E()[jsmall:jlarge+1]))
             +len(obsfwd.E()[:jsmall]))

    return jsmall,jlarge,j
