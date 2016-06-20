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
 *  $Date: 88/06/24 05:02:27 $
 *  $Revision: 1.1 $
 */




#ifndef  spOKAY
#include <stdio.h>
#include "spConfig.h"


/* Begin error macros. */
#define  spOKAY                 0
#define  spSMALL_PIVOT          1
#define  spZERO_DIAG            2
#define  spSINGULAR             3
#define  spNO_MEMORY            4
#define  spPANIC                5

#define  spFATAL                2


#if spCOMPATIBILITY
#define  NO_ERROR               spOKAY
#define  UNDER_FLOW             spOKAY
#define  OVER_FLOW              spOKAY
#define  ILL_CONDITIONED        spSMALL_PIVOT
#define  SINGULAR               spSINGULAR
#define  NO_MEMORY              spNO_MEMORY
#define  RANGE                  spPANIC

#define  FATAL                  spFATAL

#undef   spZERO_DIAG
#define  spZERO_DIAG            spSINGULAR
#endif /* spCOMPATIBILITY */


#if spCOMPATIBILITY
#define SPARSE_REAL     double
#endif


/* Begin partition keywords. */
#define spDEFAULT_PARTITION     0
#define spDIRECT_PARTITION      1
#define spINDIRECT_PARTITION    2
#define spAUTO_PARTITION        3


/* Begin Macros. */
#define  spADD_REAL_ELEMENT(element,real)       *(element) += real

#define  spADD_IMAG_ELEMENT(element,imag)       *(element+1) += imag

#define  spADD_COMPLEX_ELEMENT(element,real,imag)       \
{   *(element) += real;                                 \
    *(element+1) += imag;                               \
}

#define  spADD_REAL_QUAD(template,real)         \
{   *((template).Element1) += real;             \
    *((template).Element2) += real;             \
    *((template).Element3Negated) -= real;      \
    *((template).Element4Negated) -= real;      \
}

#define  spADD_IMAG_QUAD(template,imag)         \
{   *((template).Element1+1) += imag;           \
    *((template).Element2+1) += imag;           \
    *((template).Element3Negated+1) -= imag;    \
    *((template).Element4Negated+1) -= imag;    \
}

#define  spADD_COMPLEX_QUAD(template,real,imag) \
{   *((template).Element1) += real;             \
    *((template).Element2) += real;             \
    *((template).Element3Negated) -= real;      \
    *((template).Element4Negated) -= real;      \
    *((template).Element1+1) += imag;           \
    *((template).Element2+1) += imag;           \
    *((template).Element3Negated+1) -= imag;    \
    *((template).Element4Negated+1) -= imag;    \
}

#if spCOMPATIBILITY
#define  ADD_REAL_ELEMENT_TO_MATRIX             spADD_REAL_ELEMENT
#define  ADD_IMAG_ELEMENT_TO_MATRIX             spADD_IMAG_ELEMENT
#define  ADD_COMPLEX_ELEMENT_TO_MATRIX          spADD_COMPLEX_ELEMENT
#define  ADD_REAL_QUAD_ELEMENT_TO_MATRIX        spADD_REAL_ELEMENT
#define  ADD_IMAG_QUAD_ELEMENT_TO_MATRIX        spADD_IMAG_ELEMENT
#define  ADD_COMPLEX_QUAD_ELEMENT_TO_MATRIX     spADD_COMPLEX_ELEMENT
#endif


#if spCOMPATIBILITY
#define spTemplate TemplateStruct
#endif

/* Begin `spTemplate'. */
struct spTemplate
{
  double *Element1;
  double *Element2;
  double *Element3Negated;
  double *Element4Negated;
};


#ifdef __STDC__

/* For compilers that understand function prototypes. */

extern void spClear (char *);
extern double spCondition (char *, double, int *);
extern char *spCreate (int, int, int *);
extern void spDeleteRowAndCol (char *, int, int);
extern void spDestroy (char *);
extern int spElementCount (char *);
extern int spError (char *);
extern int spFactor (char *);
extern int spFileMatrix (char *, char *, char *, int, int, int);
extern int spFileStats (char *, char *, char *);
extern int spFillinCount (char *);
extern int spGetAdmittance (char *, int, int, struct spTemplate *);
extern double *spGetElement (char *, int, int);
// Added by hba, 23-Sep-98
extern double *spFindElement (char *, int, int);
extern char *spGetInitInfo (double *);
extern int spGetOnes (char *, int, int, int, struct spTemplate *);
extern int spGetQuad (char *, int, int, int, int, struct spTemplate *);
extern int spGetSize (char *, int);
extern int spInitialize (char *, int (*)());
extern void spInstallInitInfo (double *, char *);
extern double spLargestElement (char *);
extern void spMNA_Preorder (char *);
extern double spNorm (char *);
extern int spOrderAndFactor (char *, double[], double, double, int);
extern void spPartition (char *, int);
//hba250998 extern void spPrint (char *, int, int, int);
extern void spPrint (char *, int, int, int, FILE *);
extern double spPseudoCondition (char *);
extern double spRoundoff (char *, double);
extern void spScale (char *, double[], double[]);
extern void spSetComplex (char *);
extern void spSetReal (char *);
extern void spStripFills (char *);
extern void spWhereSingular (char *, int *, int *);

/* Functions with argument lists that are dependent on options. */

#if spCOMPLEX
extern void spDeterminant (char *, int *, double *, double *);
#else /* NOT spCOMPLEX */
extern void spDeterminant (char *, int *, double *);
#endif /* NOT spCOMPLEX */
#if spCOMPLEX && spSEPARATED_COMPLEX_VECTORS
extern int spFileVector (char *, char *, double[], double[]);
extern void spMultiply (char *, double[], double[], double[], double[]);
extern void spMultTransposed (char *, double[], double[], double[], double[]);
extern void spSolve (char *, double[], double[], double[], double[]);
extern void spSolveTransposed (char *, double[], double[], double[], double[]);
#else /* NOT (spCOMPLEX && spSEPARATED_COMPLEX_VECTORS) */
extern int spFileVector (char *, char *, double[]);
extern void spMultiply (char *, double[], double[]);
extern void spMultTransposed (char *, double[], double[]);
extern void spSolve (char *, double[], double[]);
extern void spSolveTransposed (char *, double[], double[]);
#endif /* NOT (spCOMPLEX && spSEPARATED_COMPLEX_VECTORS) */

#else /* NOT defined(__STDC__) */

/* For compilers that do not understand function prototypes. */

extern void spClear ();
extern double spCondition ();
extern char *spCreate ();
extern void spDeleteRowAndCol ();
extern void spDestroy ();
extern void spDeterminant ();
extern int spElementCount ();
extern int spError ();
extern int spFactor ();
extern int spFileMatrix ();
extern int spFileStats ();
extern int spFileVector ();
extern int spFillinCount ();
extern int spGetAdmittance ();
extern double *spGetElement ();
extern char *spGetInitInfo ();
extern int spGetOnes ();
extern int spGetQuad ();
extern int spGetSize ();
extern int spInitialize ();
extern void spInstallInitInfo ();
extern double spLargestElement ();
extern void spMNA_Preorder ();
extern void spMultiply ();
extern void spMultTransposed ();
extern double spNorm ();
extern int spOrderAndFactor ();
extern void spPartition ();
extern void spPrint ();
extern double spPseudoCondition ();
extern double spRoundoff ();
extern void spScale ();
extern void spSetComplex ();
extern void spSetReal ();
extern void spSolve ();
extern void spSolveTransposed ();
extern void spStripFills ();
extern void spWhereSingular ();
#endif /* defined(__STDC__) */

#if spCOMPATIBILITY
extern char *AllocateMatrix ();
extern double *AddElementToMatrix ();
extern void AddRealElementToMatrix ();
extern void AddImagElementToMatrix ();
extern void AddComplexElementToMatrix ();
extern void AddAdmittanceToMatrix ();
extern void AddOnesToMatrix ();
extern void AddQuadToMatrix ();
extern void AddRealQuadElementToMatrix ();
extern void AddImagQuadElementToMatrix ();
extern void AddComplexQuadElementToMatrix ();
extern void CleanMatrix ();
extern void ClearMatrix ();
extern int ClearMatrixError ();
extern void DeallocateMatrix ();
extern void DeleteRowAndColFromMatrix ();
extern void Determinant ();
extern int DecomposeMatrix ();
extern int GetMatrixSize ();
extern int MatrixElementCount ();
extern int MatrixFillinCount ();
extern void MatrixMultiply ();
extern double MatrixRoundoffError ();
extern int MatrixError ();
extern int OrderAndDecomposeMatrix ();
extern void OutputMatrixToFile ();
extern void OutputStatisticsToFile ();
extern void OutputVectorToFile ();
extern void PreorderForModifiedNodal ();
extern void PrintMatrix ();
extern void SetMatrixComplex ();
extern void SetMatrixReal ();
extern void SolveMatrix ();
extern void SolveTransposedMatrix ();
extern void ScaleMatrix ();
#endif /* spCOMPATIBILITY */

#endif /* spOKAY */

#ifdef __cplusplus
	}
#endif  
