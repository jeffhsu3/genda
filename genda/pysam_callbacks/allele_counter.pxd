import numpy as np
cimport numpy as np

#ctypedef unsigned short[:] Vector

#ctypedef np.uint16_t DTYPE_t
cdef class AlleleCounter:
    cdef char *region
    cdef int position
    cdef int phredThreshold
    cdef int A_n
    cdef int G_n
    cdef int T_n
    cdef int C_n
    #cdef Vector counts
    #cdef np.ndarray[unsigned short, ndim=1] counts
    #cdef DTYPE_t [:] counts
