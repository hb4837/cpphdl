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
 *  $Date: 88/06/18 11:13:40 $
 *  $Revision: 1.2 $
 */


#include <stdio.h>


/*
 *  If running lint, change some of the compiler options to get a more
 *  complete inspection.
 */

#ifdef lint
#undef  REAL
#undef  spCOMPLEX
#undef  EXPANDABLE
#undef  TRANSLATE
#undef  INITIALIZE
#undef  DELETE
#undef  STRIP
#undef  MODIFIED_NODAL
#undef  QUAD_ELEMENT
#undef  TRANSPOSE
#undef  SCALING
#undef  DOCUMENTATION
#undef  MULTIPLICATION
#undef  DETERMINANT
#undef  CONDITION
#undef  PSEUDOCONDITION
#undef  FORTRAN
#undef  DEBUG
#undef  spCOMPATIBILITY

#define  REAL                           YES
#define  spCOMPLEX                      YES
#define  EXPANDABLE                     YES
#define  TRANSLATE                      YES
#define  INITIALIZE                     YES
#define  DELETE                         YES
#define  STRIP                          YES
#define  MODIFIED_NODAL                 YES
#define  QUAD_ELEMENT                   YES
#define  TRANSPOSE                      YES
#define  SCALING                        YES
#define  DOCUMENTATION                  YES
#define  MULTIPLICATION                 YES
#define  DETERMINANT                    YES
#define  CONDITION                      YES
#define  PSEUDOCONDITION                YES
#define  FORTRAN                        YES
#define  DEBUG                          YES
#define  spCOMPATIBILITY                YES

#define  LINT                           YES
#else /* not lint */
#define  LINT                           NO
#endif /* not lint */


/* Begin macros. */

/* Boolean data type */
#define  BOOLEAN        int
#define  NO             0
#define  YES            1
#define  NOT            !
#define  AND            &&
#define  OR             ||

/* NULL pointer */
#ifndef  NULL
#define  NULL           0
#endif

#define  SPARSE_ID      0x772773/* Arbitrary (is Sparse on phone). */
#define  IS_SPARSE(matrix)      ((matrix) != NULL &&            \
                                 (matrix)->ID == SPARSE_ID)
#define  IS_VALID(matrix)       ((matrix) != NULL &&            \
                                 (matrix)->ID == SPARSE_ID &&   \
                                 (matrix)->Error >= spOKAY &&   \
                                 (matrix)->Error < spFATAL)
#define  IS_FACTORED(matrix)    ((matrix)->Factored && !(matrix)->NeedsOrdering)

/* Macro commands */
/* Macro functions that return the maximum or minimum independent of type. */
#define  MAX(a,b)           ((a) > (b) ? (a) : (b))
#define  MIN(a,b)           ((a) < (b) ? (a) : (b))

/* Macro function that returns the absolute value of a floating point number. */
#define  ABS(a)             ((a) < 0.0 ? -(a) : (a))

/* Macro function that returns the square of a number. */
#define  SQR(a)             ((a)*(a))

/* Macro procedure that swaps two entities. */
#define  SWAP(type, a, b)   {type swapx; swapx = a; a = b; b = swapx;}

/* Macro function that returns the approx absolute value of a complex number. */
#if spCOMPLEX
#define  ELEMENT_MAG(ptr)   (ABS((ptr)->Real) + ABS((ptr)->Imag))
#else
#define  ELEMENT_MAG(ptr)   ((ptr)->Real < 0.0 ? -(ptr)->Real : (ptr)->Real)
#endif

/* Complex assignment statements. */
#define  CMPLX_ASSIGN(to,from)  \
{   (to).Real = (from).Real;    \
    (to).Imag = (from).Imag;    \
}
#define  CMPLX_CONJ_ASSIGN(to,from)     \
{   (to).Real = (from).Real;            \
    (to).Imag = -(from).Imag;           \
}
#define  CMPLX_NEGATE_ASSIGN(to,from)   \
{   (to).Real = -(from).Real;           \
    (to).Imag = -(from).Imag;           \
}
#define  CMPLX_CONJ_NEGATE_ASSIGN(to,from)      \
{   (to).Real = -(from).Real;                   \
    (to).Imag = (from).Imag;                    \
}
#define  CMPLX_CONJ(a)  (a).Imag = -(a).Imag
#define  CMPLX_NEGATE(a)        \
{   (a).Real = -(a).Real;       \
    (a).Imag = -(a).Imag;       \
}

/* Macro that returns the approx magnitude (L-1 norm) of a complex number. */
#define  CMPLX_1_NORM(a)        (ABS((a).Real) + ABS((a).Imag))

/* Macro that returns the approx magnitude (L-infinity norm) of a complex. */
#define  CMPLX_INF_NORM(a)      (MAX (ABS((a).Real),ABS((a).Imag)))

/* Macro function that returns the magnitude (L-2 norm) of a complex number. */
#define  CMPLX_2_NORM(a)        (sqrt((a).Real*(a).Real + (a).Imag*(a).Imag))

/* Macro function that performs complex addition. */
#define  CMPLX_ADD(to,from_a,from_b)            \
{   (to).Real = (from_a).Real + (from_b).Real;  \
    (to).Imag = (from_a).Imag + (from_b).Imag;  \
}

/* Macro function that performs complex subtraction. */
#define  CMPLX_SUBT(to,from_a,from_b)           \
{   (to).Real = (from_a).Real - (from_b).Real;  \
    (to).Imag = (from_a).Imag - (from_b).Imag;  \
}

/* Macro function that is equivalent to += operator for complex numbers. */
#define  CMPLX_ADD_ASSIGN(to,from)      \
{   (to).Real += (from).Real;           \
    (to).Imag += (from).Imag;           \
}

/* Macro function that is equivalent to -= operator for complex numbers. */
#define  CMPLX_SUBT_ASSIGN(to,from)     \
{   (to).Real -= (from).Real;           \
    (to).Imag -= (from).Imag;           \
}

/* Macro function that multiplies a complex number by a scalar. */
#define  SCLR_MULT(to,sclr,cmplx)       \
{   (to).Real = (sclr) * (cmplx).Real;  \
    (to).Imag = (sclr) * (cmplx).Imag;  \
}

/* Macro function that multiply-assigns a complex number by a scalar. */
#define  SCLR_MULT_ASSIGN(to,sclr)      \
{   (to).Real *= (sclr);                \
    (to).Imag *= (sclr);                \
}

/* Macro function that multiplies two complex numbers. */
#define  CMPLX_MULT(to,from_a,from_b)           \
{   (to).Real = (from_a).Real * (from_b).Real - \
                (from_a).Imag * (from_b).Imag;  \
    (to).Imag = (from_a).Real * (from_b).Imag + \
                (from_a).Imag * (from_b).Real;  \
}

/* Macro function that implements to *= from for complex numbers. */
#define  CMPLX_MULT_ASSIGN(to,from)             \
{   RealNumber to_real_ = (to).Real;            \
    (to).Real = to_real_ * (from).Real -        \
                (to).Imag * (from).Imag;        \
    (to).Imag = to_real_ * (from).Imag +        \
                (to).Imag * (from).Real;        \
}

/* Macro function that multiplies two complex numbers, the first of which is
 * conjugated. */
#define  CMPLX_CONJ_MULT(to,from_a,from_b)      \
{   (to).Real = (from_a).Real * (from_b).Real + \
                (from_a).Imag * (from_b).Imag;  \
    (to).Imag = (from_a).Real * (from_b).Imag - \
                (from_a).Imag * (from_b).Real;  \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to = add + mult_a * mult_b */
#define  CMPLX_MULT_ADD(to,mult_a,mult_b,add)                   \
{   (to).Real = (mult_a).Real * (mult_b).Real -                 \
                (mult_a).Imag * (mult_b).Imag + (add).Real;     \
    (to).Imag = (mult_a).Real * (mult_b).Imag +                 \
                (mult_a).Imag * (mult_b).Real + (add).Imag;     \
}

/* Macro function that subtracts the product of two complex numbers from
 * another.  to = subt - mult_a * mult_b */
#define  CMPLX_MULT_SUBT(to,mult_a,mult_b,subt)                 \
{   (to).Real = (subt).Real - (mult_a).Real * (mult_b).Real +   \
                              (mult_a).Imag * (mult_b).Imag;    \
    (to).Imag = (subt).Imag - (mult_a).Real * (mult_b).Imag -   \
                              (mult_a).Imag * (mult_b).Real;    \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to = add + mult_a* * mult_b where mult_a* represents mult_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_ADD(to,mult_a,mult_b,add)              \
{   (to).Real = (mult_a).Real * (mult_b).Real +                 \
                (mult_a).Imag * (mult_b).Imag + (add).Real;     \
    (to).Imag = (mult_a).Real * (mult_b).Imag -                 \
                (mult_a).Imag * (mult_b).Real + (add).Imag;     \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to another. to += mult_a * mult_b */
#define  CMPLX_MULT_ADD_ASSIGN(to,from_a,from_b)        \
{   (to).Real += (from_a).Real * (from_b).Real -        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag += (from_a).Real * (from_b).Imag +        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then subtracts them
 * from another. */
#define  CMPLX_MULT_SUBT_ASSIGN(to,from_a,from_b)       \
{   (to).Real -= (from_a).Real * (from_b).Real -        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag -= (from_a).Real * (from_b).Imag +        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then adds them
 * to the destination. to += from_a* * from_b where from_a* represents from_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_ADD_ASSIGN(to,from_a,from_b)   \
{   (to).Real += (from_a).Real * (from_b).Real +        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag += (from_a).Real * (from_b).Imag -        \
                 (from_a).Imag * (from_b).Real;         \
}

/* Macro function that multiplies two complex numbers and then subtracts them
 * from the destination. to -= from_a* * from_b where from_a* represents from_a
 * conjugate. */
#define  CMPLX_CONJ_MULT_SUBT_ASSIGN(to,from_a,from_b)  \
{   (to).Real -= (from_a).Real * (from_b).Real +        \
                 (from_a).Imag * (from_b).Imag;         \
    (to).Imag -= (from_a).Real * (from_b).Imag -        \
                 (from_a).Imag * (from_b).Real;         \
}

/*
 * Macro functions that provide complex division.
 */

/* Complex division:  to = num / den */
#define CMPLX_DIV(to,num,den)                                           \
{   RealNumber  r_, s_;                                                 \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        s_ = (den).Real + r_*(den).Imag;                                \
        (to).Real = ((num).Real + r_*(num).Imag)/s_;                    \
        (to).Imag = ((num).Imag - r_*(num).Real)/s_;                    \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        s_ = (den).Imag + r_*(den).Real;                                \
        (to).Real = (r_*(num).Real + (num).Imag)/s_;                    \
        (to).Imag = (r_*(num).Imag - (num).Real)/s_;                    \
    }                                                                   \
}

/* Complex division and assignment:  num /= den */
#define CMPLX_DIV_ASSIGN(num,den)                                       \
{   RealNumber  r_, s_, t_;                                             \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        s_ = (den).Real + r_*(den).Imag;                                \
        t_ = ((num).Real + r_*(num).Imag)/s_;                           \
        (num).Imag = ((num).Imag - r_*(num).Real)/s_;                   \
        (num).Real = t_;                                                \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        s_ = (den).Imag + r_*(den).Real;                                \
        t_ = (r_*(num).Real + (num).Imag)/s_;                           \
        (num).Imag = (r_*(num).Imag - (num).Real)/s_;                   \
        (num).Real = t_;                                                \
    }                                                                   \
}

/* Complex reciprocation:  to = 1.0 / den */
#define CMPLX_RECIPROCAL(to,den)                                        \
{   RealNumber  r_;                                                     \
    if (((den).Real >= (den).Imag AND (den).Real > -(den).Imag) OR      \
        ((den).Real < (den).Imag AND (den).Real <= -(den).Imag))        \
    {   r_ = (den).Imag / (den).Real;                                   \
        (to).Imag = -r_*((to).Real = 1.0/((den).Real + r_*(den).Imag)); \
    }                                                                   \
    else                                                                \
    {   r_ = (den).Real / (den).Imag;                                   \
        (to).Real = -r_*((to).Imag = -1.0/((den).Imag + r_*(den).Real));\
    }                                                                   \
}


#if DEBUG
#define ASSERT(condition) if (NOT(condition)) ABORT()
#else
#define ASSERT(condition)
#endif

#if DEBUG
#define  ABORT()                                                        \
{   (void)fflush(stdout);                                               \
    (void)fprintf(stderr, "sparse: panic in file `%s' at line %d.\n",   \
            __FILE__, __LINE__);                                        \
    (void)fflush(stderr);                                               \
    abort();                                                            \
}
#else
#define  ABORT()
#endif


#if spCOMPLEX AND spSEPARATED_COMPLEX_VECTORS
#define IMAG_VECTORS    , iRHS, iSolution
#define IMAG_RHS        , iRHS
#else
#define IMAG_VECTORS
#define IMAG_RHS
#endif


//extern char *malloc (), *calloc (), *realloc ();
#ifdef ultrix
extern void free ();
extern void abort ();
#elif sun
extern void abort();
//#else
//extern free ();
//extern abort ();
#endif

#define ALLOC(type,number)  ((type *)malloc((unsigned)(sizeof(type)*(number))))
#define REALLOC(ptr,type,number)  \
           ptr = (type *)realloc((char *)ptr,(unsigned)(sizeof(type)*(number)))
#define FREE(ptr) { if ((ptr) != NULL) free((char *)(ptr)); (ptr) = NULL; }


/* Calloc that properly handles allocating a cleared vector. */
#define CALLOC(ptr,type,number)                         \
{   int i; ptr = ALLOC(type, number);                   \
    if (ptr != (type *)NULL)                            \
        for(i=(number)-1;i>=0; i--) ptr[i] = (type) 0;  \
}


/* Begin `RealNumber'. */
typedef double RealNumber, *RealVector;


/* Begin `ComplexNumber'. */
typedef struct
{
  RealNumber Real;
  RealNumber Imag;
} ComplexNumber, *ComplexVector;


/* Begin `MatrixElement'. */
struct MatrixElement
{
  RealNumber Real;
#if spCOMPLEX
  RealNumber Imag;
#endif
  int Row;
  int Col;
  struct MatrixElement *NextInRow;
  struct MatrixElement *NextInCol;
#if INITIALIZE
  char *pInitInfo;
#endif
};

typedef struct MatrixElement *ElementPtr;
typedef ElementPtr *ArrayOfElementPtrs;


struct AllocationRecord
{
  char *AllocatedPtr;
  struct AllocationRecord *NextRecord;
};

typedef struct AllocationRecord *AllocationListPtr;


struct FillinListNodeStruct
{
  ElementPtr pFillinList;
  int NumberOfFillinsInList;
  struct FillinListNodeStruct *Next;
};


struct MatrixFrame
{
  RealNumber AbsThreshold;
  int AllocatedSize;
  int AllocatedExtSize;
  BOOLEAN Complex;
  int CurrentSize;
  ArrayOfElementPtrs Diag;
  BOOLEAN *DoCmplxDirect;
  BOOLEAN *DoRealDirect;
  int Elements;
  int Error;
  int ExtSize;
  int *ExtToIntColMap;
  int *ExtToIntRowMap;
  BOOLEAN Factored;
  int Fillins;
  ArrayOfElementPtrs FirstInCol;
  ArrayOfElementPtrs FirstInRow;
  unsigned long ID;
  RealVector Intermediate;
  BOOLEAN InternalVectorsAllocated;
  int *IntToExtColMap;
  int *IntToExtRowMap;
  int *MarkowitzRow;
  int *MarkowitzCol;
  long *MarkowitzProd;
  int MaxRowCountInLowerTri;
  BOOLEAN NeedsOrdering;
  BOOLEAN NumberOfInterchangesIsOdd;
  BOOLEAN Partitioned;
  int PivotsOriginalCol;
  int PivotsOriginalRow;
  char PivotSelectionMethod;
  BOOLEAN PreviousMatrixWasComplex;
  RealNumber RelThreshold;
  BOOLEAN Reordered;
  BOOLEAN RowsLinked;
  int SingularCol;
  int SingularRow;
  int Singletons;
  int Size;
  struct MatrixElement TrashCan;

  AllocationListPtr TopOfAllocationList;
  int RecordsRemaining;
  ElementPtr NextAvailElement;
  int ElementsRemaining;
  ElementPtr NextAvailFillin;
  int FillinsRemaining;
  struct FillinListNodeStruct *FirstFillinListNode;
  struct FillinListNodeStruct *LastFillinListNode;
};
typedef struct MatrixFrame *MatrixPtr;

#ifdef __cplusplus
	}
#endif  
