#ifdef __cplusplus
	extern "C" {
#endif  
/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 *
 *  $Date: 88/06/24 05:00:58 $
 *  $Revision: 1.3 $
 */


#ifndef spCONFIG_DEFS
#define spCONFIG_DEFS




#ifdef spINSIDE_SPARSE

/* Begin options. */
#define  REAL                           YES
#define  EXPANDABLE                     YES
#define  TRANSLATE                      YES
#define  INITIALIZE                     YES
#define  DIAGONAL_PIVOTING              YES
#define  ARRAY_OFFSET                   NO
#define  MODIFIED_MARKOWITZ             NO
#define  DELETE                         YES
#define  STRIP                          YES
#define  MODIFIED_NODAL                 YES
#define  QUAD_ELEMENT                   YES
#define  TRANSPOSE                      YES
#define  SCALING                        YES
#define  DOCUMENTATION                  YES
#define  MULTIPLICATION                 YES
#define  DETERMINANT                    YES
#define  STABILITY                      YES
#define  CONDITION                      YES
#define  PSEUDOCONDITION                YES
#define  FORTRAN                        NO
//#define  DEBUG                          YES

/*
 *  The following options affect Sparse exports and so are exported as a
 *  side effect.  For this reason they use the `sp' prefix.  The boolean
 *  constants YES an NO are not defined in spMatrix.h to avoid conflicts
 *  with user code, so use 0 for NO and 1 for YES.
 */
#endif /* spINSIDE_SPARSE */
#define  spCOMPLEX                      1
#define  spSEPARATED_COMPLEX_VECTORS    0
#define  spCOMPATIBILITY                0
#ifdef spINSIDE_SPARSE


/* Begin constants. */
#define  DEFAULT_THRESHOLD              1.0e-3
#define  DIAG_PIVOTING_AS_DEFAULT       YES
#define  SPACE_FOR_ELEMENTS             6
#define  SPACE_FOR_FILL_INS             4
#define  ELEMENTS_PER_ALLOCATION        31
#define  MINIMUM_ALLOCATED_SIZE         6
#define  EXPANSION_FACTOR               1.5
#define  MAX_MARKOWITZ_TIES             100
#define  TIES_MULTIPLIER                5
#define  DEFAULT_PARTITION              spAUTO_PARTITION


/*  Begin printer constants. */
#define  PRINTER_WIDTH  80


/* Begin machine constants. */

#ifdef notdef			/* __STDC__ */
/*
 * This code is currently deleted because most ANSI standard C compilers
 * do not provide the standard header files yet.
 */
#   include <limits.h>
#   include <float.h>
#   define  MACHINE_RESOLUTION      DBL_EPSILON
#   define  LARGEST_REAL            DBL_MAX
#   define  SMALLEST_REAL           DBL_MIN
#   define  LARGEST_SHORT_INTEGER   SHRT_MAX
#   define  LARGEST_LONG_INTEGER    LONG_MAX
#else /* NOT defined(__STDC__) */

/* VAX machine constants */
#ifdef vax
#   define  MACHINE_RESOLUTION      6.93889e-18
#   define  LARGEST_REAL            1.70141e+38
#   define  SMALLEST_REAL           2.938743e-39
#   define  LARGEST_SHORT_INTEGER   32766
#   define  LARGEST_LONG_INTEGER    2147483646
#endif

/* hp9000 machine constants */
#ifdef hpux
/* These values are correct for hp9000/300.  Should be correct for others. */
#   define  MACHINE_RESOLUTION      8.9e-15
#   define  LARGEST_REAL            1.79769313486231e+308
#   define  SMALLEST_REAL           2.22507385850721e-308
#   define  LARGEST_SHORT_INTEGER   32766
#   define  LARGEST_LONG_INTEGER    2147483646
#endif

/* Sun machine constants */
#ifdef sun
/* These values are rumored to be the correct values. */
#   define  MACHINE_RESOLUTION      8.9e-15
#   define  LARGEST_REAL            1.79769313486231e+308
#   define  SMALLEST_REAL           2.22507385850721e-308
#   define  LARGEST_SHORT_INTEGER   32766
#   define  LARGEST_LONG_INTEGER    2147483646
#endif

#ifdef linux
#   define  MACHINE_RESOLUTION      8.9e-15
#   define  LARGEST_REAL            1.79769313486231e+308
#   define  SMALLEST_REAL           2.22507385850721e-308
#   define  LARGEST_SHORT_INTEGER   32766
#   define  LARGEST_LONG_INTEGER    2147483646
#endif


#ifdef _WIN32
#   define  MACHINE_RESOLUTION      8.9e-15
#   define  LARGEST_REAL            1.79769313486231e+308
#   define  SMALLEST_REAL           2.22507385850721e-308
#   define  LARGEST_SHORT_INTEGER   32766
#   define  LARGEST_LONG_INTEGER    2147483646
#endif
#endif /* NOT defined(__STDC__) */


/* Begin annotation definitions. */
#define  ANNOTATE               SPNONE

#define  SPNONE                   0
#define  ON_STRANGE_BEHAVIOR    1
#define  FULL                   2

#endif /* spINSIDE_SPARSE */
#endif /* spCONFIG_DEFS */

#ifdef __cplusplus
	}
#endif  
