set(NSOURCES
        denseMatrix.cc
        nnls.cc
        sparseMatrix.cc
        )

set(NHEADERS
        denseMatrix.h
        matrix.h
        nnls.h
        sparseMatrix.h
        vector.h
        )

add_library(nnls_lib STATIC ${NSOURCES} ${NHEADERS})

