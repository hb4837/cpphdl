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
 */

#ifndef lint
static char copyright[] =
    "Sparse1.3: Copyright (c) 1985,86,87,88 by Kenneth S. Kundert";
static char RCSid[] =
    "@(#)$Header: spUtils.c,v 1.3 88/06/24 05:03:37 kundert Exp $";
#endif


#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"


static int CountTwins(MatrixPtr, int, ElementPtr *, ElementPtr *);
static SwapCols(MatrixPtr, ElementPtr, ElementPtr);
static void ScaleComplexMatrix(MatrixPtr, RealVector, RealVector);
static void ComplexMatrixMultiply(MatrixPtr Matrix, RealVector RHS, RealVector Solution IMAG_VECTORS);
static void ComplexTransposedMatrixMultiply(MatrixPtr, RealVector, RealVector Solution IMAG_VECTORS);
static double ComplexCondition(MatrixPtr, double, int *);

#if MODIFIED_NODAL
void spMNA_Preorder( eMatrix )
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  int  J, Size;
ElementPtr  pTwin1, pTwin2;
int  Twins, StartAt = 1;
BOOLEAN  Swapped, AnotherPassNeeded;

/* Begin `spMNA_Preorder'. */
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored );

    if (Matrix->RowsLinked) return;
    Size = Matrix->Size;
    Matrix->Reordered = YES;

    do
    {   AnotherPassNeeded = Swapped = NO;

/* Search for zero diagonals with lone twins. */
        for (J = StartAt; J <= Size; J++)
        {   if (Matrix->Diag[J] == NULL)
            {   Twins = CountTwins( Matrix, J, &pTwin1, &pTwin2 );
                if (Twins == 1)
                {   /* Lone twins found, swap rows. */
                    SwapCols( Matrix, pTwin1, pTwin2 );
                    Swapped = YES;
                }
                else if ((Twins > 1) AND NOT AnotherPassNeeded)
                {   AnotherPassNeeded = YES;
                    StartAt = J;
                }
            }
        }

/* All lone twins are gone, look for zero diagonals with multiple twins. */
        if (AnotherPassNeeded)
        {   for (J = StartAt; NOT Swapped AND (J <= Size); J++)
            {   if (Matrix->Diag[J] == NULL)
                {   Twins = CountTwins( Matrix, J, &pTwin1, &pTwin2 );
                    SwapCols( Matrix, pTwin1, pTwin2 );
                    Swapped = YES;
                }
            }
        }
    } while (AnotherPassNeeded);
    return;
}


static int CountTwins( Matrix, Col, ppTwin1, ppTwin2 )
MatrixPtr Matrix;
int Col;
ElementPtr *ppTwin1, *ppTwin2;
{
int Row, Twins = 0;
ElementPtr pTwin1, pTwin2;

/* Begin `CountTwins'. */

    pTwin1 = Matrix->FirstInCol[Col];
    while (pTwin1 != NULL)
    {   if (ABS(pTwin1->Real) == 1.0)
        {   Row = pTwin1->Row;
            pTwin2 = Matrix->FirstInCol[Row];
            while ((pTwin2 != NULL) AND (pTwin2->Row != Col))
                pTwin2 = pTwin2->NextInCol;
            if ((pTwin2 != NULL) AND (ABS(pTwin2->Real) == 1.0))
            {   /* Found symmetric twins. */
                if (++Twins >= 2) return Twins;
                (*ppTwin1 = pTwin1)->Col = Col;
                (*ppTwin2 = pTwin2)->Col = Row;
            }
        }
        pTwin1 = pTwin1->NextInCol;
    }
    return Twins;
}




/*
 *  SWAP COLUMNS
 *
 *  This function swaps two columns and is applicable before the rows are
 *  linked.
 */

static
SwapCols( Matrix, pTwin1, pTwin2 )

MatrixPtr Matrix;
ElementPtr pTwin1, pTwin2;
{
int Col1 = pTwin1->Col, Col2 = pTwin2->Col;

/* Begin `SwapCols'. */

    SWAP (ElementPtr, Matrix->FirstInCol[Col1], Matrix->FirstInCol[Col2]);
    SWAP (int, Matrix->IntToExtColMap[Col1], Matrix->IntToExtColMap[Col2]);
#if TRANSLATE
    Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col2]]=Col2;
    Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col1]]=Col1;
#endif

    Matrix->Diag[Col1] = pTwin2;
    Matrix->Diag[Col2] = pTwin1;
    Matrix->NumberOfInterchangesIsOdd = NOT Matrix->NumberOfInterchangesIsOdd;
    return;
}
#endif /* MODIFIED_NODAL */


#if SCALING
void spScale( eMatrix, RHS_ScaleFactors, SolutionScaleFactors )
char *eMatrix;
register  RealVector  RHS_ScaleFactors, SolutionScaleFactors;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr  pElement;
register int  I, lSize, *pExtOrder;
double  ScaleFactor;
void ScaleComplexMatrix();

/* Begin `spScale'. */
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored );
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

#if spCOMPLEX
    if (Matrix->Complex)
    {   ScaleComplexMatrix( Matrix, RHS_ScaleFactors, SolutionScaleFactors );
        return;
    }
#endif

#if REAL
    lSize = Matrix->Size;

/* Correct pointers to arrays for ARRAY_OFFSET */
#if NOT ARRAY_OFFSET
    --RHS_ScaleFactors;
    --SolutionScaleFactors;
#endif

/* Scale Rows */
    pExtOrder = &Matrix->IntToExtRowMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = RHS_ScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInRow[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement = pElement->NextInRow;
            }
        }
    }

/* Scale Columns */
    pExtOrder = &Matrix->IntToExtColMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = SolutionScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInCol[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement = pElement->NextInCol;
            }
        }
    }
    return;

#endif /* REAL */
}
#endif /* SCALING */


#if spCOMPLEX AND SCALING
static void ScaleComplexMatrix( Matrix, RHS_ScaleFactors, SolutionScaleFactors )
MatrixPtr  Matrix;
register  RealVector  RHS_ScaleFactors, SolutionScaleFactors;
{
register ElementPtr  pElement;
register int  I, lSize, *pExtOrder;
double  ScaleFactor;

/* Begin `ScaleComplexMatrix'. */
    lSize = Matrix->Size;

/* Correct pointers to arrays for ARRAY_OFFSET */
#if NOT ARRAY_OFFSET
    --RHS_ScaleFactors;
    --SolutionScaleFactors;
#endif

/* Scale Rows */
    pExtOrder = &Matrix->IntToExtRowMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = RHS_ScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInRow[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement->Imag *= ScaleFactor;
                pElement = pElement->NextInRow;
            }
        }
    }

/* Scale Columns */
    pExtOrder = &Matrix->IntToExtColMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = SolutionScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInCol[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement->Imag *= ScaleFactor;
                pElement = pElement->NextInCol;
            }
        }
    }
    return;
}
#endif /* SCALING AND spCOMPLEX */


#if MULTIPLICATION
void spMultiply( eMatrix, RHS, Solution IMAG_VECTORS )
char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  RealVector  Vector;
register  double  Sum;
register  int  I, *pExtOrder;
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
extern void ComplexMatrixMultiply();

/* Begin `spMultiply'. */
    ASSERT( IS_SPARSE( Matrix ) AND NOT Matrix->Factored );
    if (NOT Matrix->RowsLinked) spcLinkRows(Matrix);

#if spCOMPLEX
    if (Matrix->Complex)
    {   ComplexMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS );
        return;
    }
#endif

#if REAL
#if NOT ARRAY_OFFSET
/* Correct array pointers for ARRAY_OFFSET. */
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = Solution[*(pExtOrder--)];

    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInRow[I];
        Sum = 0.0;

        while (pElement != NULL)
        {   Sum += pElement->Real * Vector[pElement->Col];
            pElement = pElement->NextInRow;
        }
        RHS[*pExtOrder--] = Sum;
    }
    return;
#endif /* REAL */
}
#endif /* MULTIPLICATION */


#if spCOMPLEX AND MULTIPLICATION
static void ComplexMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS )
MatrixPtr  Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Vector;
register  ComplexNumber  Sum;
register  int  I, *pExtOrder;

/* Begin `ComplexMatrixMultiply'. */

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;              --iRHS;
    --Solution;         --iSolution;
#else
    RHS -= 2;           Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = (ComplexVector)Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Matrix->Size; I > 0; I--)
    {   Vector[I].Real = Solution[*pExtOrder];
        Vector[I].Imag = iSolution[*(pExtOrder--)];
    }
#else
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = ((ComplexVector)Solution)[*(pExtOrder--)];
#endif

    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInRow[I];
        Sum.Real = Sum.Imag = 0.0;

        while (pElement != NULL)
        {   /* Cmplx expression : Sum += Element * Vector[Col] */
            CMPLX_MULT_ADD_ASSIGN( Sum, *pElement, Vector[pElement->Col] );
            pElement = pElement->NextInRow;
        }

#if spSEPARATED_COMPLEX_VECTORS
        RHS[*pExtOrder] = Sum.Real;
        iRHS[*pExtOrder--] = Sum.Imag;
#else
        ((ComplexVector)RHS)[*pExtOrder--] = Sum;
#endif
    }
    return;
}
#endif /* spCOMPLEX AND MULTIPLICATION */


#if MULTIPLICATION AND TRANSPOSE
void spMultTransposed( eMatrix, RHS, Solution IMAG_VECTORS )
char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  RealVector  Vector;
register  double  Sum;
register  int  I, *pExtOrder;
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
extern void ComplexTransposedMatrixMultiply();

/* Begin `spMultTransposed'. */
    ASSERT( IS_SPARSE( Matrix ) AND NOT Matrix->Factored );

#if spCOMPLEX
    if (Matrix->Complex)
    {   ComplexTransposedMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS );
        return;
    }
#endif

#if REAL
#if NOT ARRAY_OFFSET
/* Correct array pointers for ARRAY_OFFSET. */
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = Solution[*(pExtOrder--)];

    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInCol[I];
        Sum = 0.0;

        while (pElement != NULL)
        {   Sum += pElement->Real * Vector[pElement->Row];
            pElement = pElement->NextInCol;
        }
        RHS[*pExtOrder--] = Sum;
    }
    return;
#endif /* REAL */
}
#endif /* MULTIPLICATION AND TRANSPOSE */


#if spCOMPLEX AND MULTIPLICATION AND TRANSPOSE
static void ComplexTransposedMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS )
MatrixPtr  Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Vector;
register  ComplexNumber  Sum;
register  int  I, *pExtOrder;

/* Begin `ComplexMatrixMultiply'. */

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;              --iRHS;
    --Solution;         --iSolution;
#else
    RHS -= 2;           Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = (ComplexVector)Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Matrix->Size; I > 0; I--)
    {   Vector[I].Real = Solution[*pExtOrder];
        Vector[I].Imag = iSolution[*(pExtOrder--)];
    }
#else
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = ((ComplexVector)Solution)[*(pExtOrder--)];
#endif

    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInCol[I];
        Sum.Real = Sum.Imag = 0.0;

        while (pElement != NULL)
        {   /* Cmplx expression : Sum += Element * Vector[Row] */
            CMPLX_MULT_ADD_ASSIGN( Sum, *pElement, Vector[pElement->Row] );
            pElement = pElement->NextInCol;
        }

#if spSEPARATED_COMPLEX_VECTORS
        RHS[*pExtOrder] = Sum.Real;
        iRHS[*pExtOrder--] = Sum.Imag;
#else
        ((ComplexVector)RHS)[*pExtOrder--] = Sum;
#endif
    }
    return;
}
#endif /* spCOMPLEX AND MULTIPLICATION AND TRANSPOSE */


#if DETERMINANT
#if spCOMPLEX 
void spDeterminant( eMatrix, pExponent, pDeterminant, piDeterminant )
double *piDeterminant;
#else
void spDeterminant( eMatrix, pExponent, pDeterminant )
#endif
char *eMatrix;
register  double *pDeterminant;
int  *pExponent;
{
register MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I, Size;
double Norm, nr, ni;
ComplexNumber Pivot, cDeterminant;

#define  NORM(a)     (nr = ABS((a).Real), ni = ABS((a).Imag), MAX (nr,ni))

/* Begin `spDeterminant'. */
    ASSERT( IS_SPARSE( Matrix ) AND IS_FACTORED(Matrix) );
    *pExponent = 0;

    if (Matrix->Error == spSINGULAR)
    {   *pDeterminant = 0.0;
#if spCOMPLEX
        if (Matrix->Complex) *piDeterminant = 0.0;
#endif
        return;
    }

    Size = Matrix->Size;
    I = 0;

#if spCOMPLEX
    if (Matrix->Complex)        /* Complex Case. */
    {   cDeterminant.Real = 1.0;
        cDeterminant.Imag = 0.0;

        while (++I <= Size)
        {   CMPLX_RECIPROCAL( Pivot, *Matrix->Diag[I] );
            CMPLX_MULT_ASSIGN( cDeterminant, Pivot );

/* Scale Determinant. */
            Norm = NORM( cDeterminant );
            if (Norm != 0.0)
            {   while (Norm >= 1.0e12)
                {   cDeterminant.Real *= 1.0e-12;
                    cDeterminant.Imag *= 1.0e-12;
                    *pExponent += 12;
                    Norm = NORM( cDeterminant );
                }
                while (Norm < 1.0e-12)
                {   cDeterminant.Real *= 1.0e12;
                    cDeterminant.Imag *= 1.0e12;
                    *pExponent -= 12;
                    Norm = NORM( cDeterminant );
                }
            }
        }

/* Scale Determinant again, this time to be between 1.0 <= x < 10.0. */
        Norm = NORM( cDeterminant );
        if (Norm != 0.0)
        {   while (Norm >= 10.0)
            {   cDeterminant.Real *= 0.1;
                cDeterminant.Imag *= 0.1;
                (*pExponent)++;
                Norm = NORM( cDeterminant );
            }
            while (Norm < 1.0)
            {   cDeterminant.Real *= 10.0;
                cDeterminant.Imag *= 10.0;
                (*pExponent)--;
                Norm = NORM( cDeterminant );
            }
        }
        if (Matrix->NumberOfInterchangesIsOdd)
            CMPLX_NEGATE( cDeterminant );
        
        *pDeterminant = cDeterminant.Real;
        *piDeterminant = cDeterminant.Imag;
    }
#endif /* spCOMPLEX */
#if REAL AND spCOMPLEX
    else
#endif
#if REAL
    {   /* Real Case. */
        *pDeterminant = 1.0;

        while (++I <= Size)
        {   *pDeterminant /= Matrix->Diag[I]->Real;

/* Scale Determinant. */
            if (*pDeterminant != 0.0)
            {   while (ABS(*pDeterminant) >= 1.0e12)
                {   *pDeterminant *= 1.0e-12;
                    *pExponent += 12;
                }
                while (ABS(*pDeterminant) < 1.0e-12)
                {   *pDeterminant *= 1.0e12;
                    *pExponent -= 12;
                }
            }
        }

/* Scale Determinant again, this time to be between 1.0 <= x < 10.0. */
        if (*pDeterminant != 0.0)
        {   while (ABS(*pDeterminant) >= 10.0)
            {   *pDeterminant *= 0.1;
                (*pExponent)++;
            }
            while (ABS(*pDeterminant) < 1.0)
            {   *pDeterminant *= 10.0;
                (*pExponent)--;
            }
        }
        if (Matrix->NumberOfInterchangesIsOdd)
            *pDeterminant = -*pDeterminant;
    }
#endif /* REAL */
}
#endif /* DETERMINANT */


#if STRIP
void spStripFills( eMatrix )
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
struct FillinListNodeStruct  *pListNode;

/* Begin `spStripFills'. */
    ASSERT( IS_SPARSE( Matrix ) );
    if (Matrix->Fillins == 0) return;
    Matrix->NeedsOrdering = YES;
    Matrix->Elements -= Matrix->Fillins;
    Matrix->Fillins = 0;

/* Mark the fill-ins. */
    {   register  ElementPtr  pFillin, pLastFillin;

        pListNode = Matrix->LastFillinListNode = Matrix->FirstFillinListNode;
        Matrix->FillinsRemaining = pListNode->NumberOfFillinsInList;
        Matrix->NextAvailFillin = pListNode->pFillinList;

        while (pListNode != NULL)
        {   pFillin = pListNode->pFillinList;
            pLastFillin = &(pFillin[ pListNode->NumberOfFillinsInList - 1 ]);
            while (pFillin <= pLastFillin)
                (pFillin++)->Row = 0;
            pListNode = pListNode->Next;
        }
    }

/* Unlink fill-ins by searching for elements marked with Row = 0. */
    {   register  ElementPtr pElement, *ppElement;
        register  int  I, Size = Matrix->Size;

/* Unlink fill-ins in all columns. */
        for (I = 1; I <= Size; I++)
        {   ppElement = &(Matrix->FirstInCol[I]);
            while ((pElement = *ppElement) != NULL)
            {   if (pElement->Row == 0)
                {   *ppElement = pElement->NextInCol;  /* Unlink fill-in. */
                    if (Matrix->Diag[pElement->Col] == pElement)
                        Matrix->Diag[pElement->Col] = NULL;
                }
                else
                    ppElement = &pElement->NextInCol;  /* Skip element. */
            }
        }

/* Unlink fill-ins in all rows. */
        for (I = 1; I <= Size; I++)
        {   ppElement = &(Matrix->FirstInRow[I]);
            while ((pElement = *ppElement) != NULL)
            {   if (pElement->Row == 0)
                    *ppElement = pElement->NextInRow;  /* Unlink fill-in. */
                else
                    ppElement = &pElement->NextInRow;  /* Skip element. */
            }
        }
    }
    return;
}
#endif


#if TRANSLATE AND DELETE
void spDeleteRowAndCol( eMatrix, Row, Col )
char *eMatrix;
int  Row, Col;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement, *ppElement, pLastElement;
int  Size, ExtRow, ExtCol;
ElementPtr  spcFindElementInCol();

/* Begin `spDeleteRowAndCol'. */

    ASSERT( IS_SPARSE(Matrix) AND Row > 0 AND Col > 0 );
    ASSERT( Row <= Matrix->ExtSize AND Col <= Matrix->ExtSize );

    Size = Matrix->Size;
    ExtRow = Row;
    ExtCol = Col;
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

    Row = Matrix->ExtToIntRowMap[Row];
    Col = Matrix->ExtToIntColMap[Col];
    ASSERT( Row > 0 AND Col > 0 );

/* Move Row so that it is the last row in the matrix. */
    if (Row != Size) spcRowExchange( Matrix, Row, Size );

/* Move Col so that it is the last column in the matrix. */
    if (Col != Size) spcColExchange( Matrix, Col, Size );

/* Correct Diag pointers. */
    if (Row == Col)
        SWAP( ElementPtr, Matrix->Diag[Row], Matrix->Diag[Size] )
    else
    {   Matrix->Diag[Row] = spcFindElementInCol( Matrix, Matrix->FirstInCol+Row,
                                                 Row, Row, NO );
        Matrix->Diag[Col] = spcFindElementInCol( Matrix, Matrix->FirstInCol+Col,
                                                 Col, Col, NO );
    }

/*
 * Delete last row and column of the matrix.
 */
/* Break the column links to every element in the last row. */
    pLastElement = Matrix->FirstInRow[ Size ];
    while (pLastElement != NULL)
    {   ppElement = &(Matrix->FirstInCol[ pLastElement->Col ]);
        while ((pElement = *ppElement) != NULL)
        {   if (pElement == pLastElement)
                *ppElement = NULL;  /* Unlink last element in column. */
            else
                ppElement = &pElement->NextInCol;  /* Skip element. */
        }
        pLastElement = pLastElement->NextInRow;
    }

/* Break the row links to every element in the last column. */
    pLastElement = Matrix->FirstInCol[ Size ];
    while (pLastElement != NULL)
    {   ppElement = &(Matrix->FirstInRow[ pLastElement->Row ]);
        while ((pElement = *ppElement) != NULL)
        {   if (pElement == pLastElement)
                *ppElement = NULL;  /* Unlink last element in row. */
            else
                ppElement = &pElement->NextInRow;  /* Skip element. */
        }
        pLastElement = pLastElement->NextInCol;
    }

/* Clean up some details. */
    Matrix->Size = Size - 1;
    Matrix->Diag[Size] = NULL;
    Matrix->FirstInRow[Size] = NULL;
    Matrix->FirstInCol[Size] = NULL;
    Matrix->CurrentSize--;
    Matrix->ExtToIntRowMap[ExtRow] = -1;
    Matrix->ExtToIntColMap[ExtCol] = -1;
    Matrix->NeedsOrdering = YES;

    return;
}
#endif


#if PSEUDOCONDITION
double spPseudoCondition( eMatrix )
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I;
register ArrayOfElementPtrs Diag;
double MaxPivot, MinPivot, Mag;

/* Begin `spPseudoCondition'. */

    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );
    if (Matrix->Error == spSINGULAR OR Matrix->Error == spZERO_DIAG)
        return 0.0;

    Diag = Matrix->Diag;
    MaxPivot = MinPivot = ELEMENT_MAG( Diag[1] );
    for (I = 2; I <= Matrix->Size; I++)
    {   Mag = ELEMENT_MAG( Diag[I] );
        if (Mag > MaxPivot)
            MaxPivot = Mag;
        else if (Mag < MinPivot)
            MinPivot = Mag;
    }
    ASSERT( MaxPivot > 0.0);
    return MaxPivot / MinPivot;
}
#endif


#if CONDITION
double spCondition( eMatrix, NormOfMatrix, pError )
char *eMatrix;
double NormOfMatrix;
int *pError;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register RealVector T, Tm;
register int I, K, Row;
ElementPtr pPivot;
int Size;
double E, Em, Wp, Wm, ASp, ASm, ASw, ASy, ASv, ASz, MaxY, ScaleFactor;
double Linpack, OLeary, InvNormOfInverse, ComplexCondition();
#define SLACK   1e4

/* Begin `spCondition'. */

    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );
    *pError = Matrix->Error;
    if (Matrix->Error >= spFATAL) return 0.0;
    if (NormOfMatrix == 0.0)
    {   *pError = spSINGULAR;
        return 0.0;
    }

#if spCOMPLEX
    if (Matrix->Complex)
        return ComplexCondition( Matrix, NormOfMatrix, pError );
#endif

#if REAL
    Size = Matrix->Size;
    T = Matrix->Intermediate;
#if spCOMPLEX
    Tm = Matrix->Intermediate + Size;
#else
    Tm = ALLOC( double, Size+1 );
    if (Tm == NULL)
    {   *pError = spNO_MEMORY;
        return 0.0;
    }
#endif
    for (I = Size; I > 0; I--) T[I] = 0.0;

/*
 * Part 1.  Ay = e.
 * Solve Ay = LUy = e where e consists of +1 and -1 terms with the sign
 * chosen to maximize the size of w in Lw = e.  Since the terms in w can
 * get very large, scaling is used to avoid overflow.
 */

/* Forward elimination. Solves Lw = e while choosing e. */
    E = 1.0;
    for (I = 1; I <= Size; I++)
    {   pPivot = Matrix->Diag[I];
        if (T[I] < 0.0) Em = -E; else Em = E;
        Wm = (Em + T[I]) * pPivot->Real;
        if (ABS(Wm) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(Wm) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            E *= ScaleFactor;
            Em *= ScaleFactor;
            Wm = (Em + T[I]) * pPivot->Real;
        }
        Wp = (T[I] - Em) * pPivot->Real;
        ASp = ABS(T[I] - Em);
        ASm = ABS(Em + T[I]);

/* Update T for both values of W, minus value is placed in Tm. */
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   Row = pElement->Row;
            Tm[Row] = T[Row] - (Wm * pElement->Real);
            T[Row] -= (Wp * pElement->Real);
            ASp += ABS(T[Row]);
            ASm += ABS(Tm[Row]);
            pElement = pElement->NextInCol;
        }

/* If minus value causes more growth, overwrite T with its values. */
        if (ASm > ASp)
        {   T[I] = Wm;
            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {   T[pElement->Row] = Tm[pElement->Row];
                pElement = pElement->NextInCol;
            }
        }
        else T[I] = Wp;
    }

/* Compute 1-norm of T, which now contains w, and scale ||T|| to 1/SLACK. */
    for (ASw = 0.0, I = Size; I > 0; I--) ASw += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASw);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        E *= ScaleFactor;
    }

/* Backward Substitution. Solves Uy = w.*/
    for (I = Size; I >= 1; I--)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   T[I] -= pElement->Real * T[pElement->Col];
            pElement = pElement->NextInRow;
        }
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            E *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains y, and scale ||T|| to 1/SLACK. */
    for (ASy = 0.0, I = Size; I > 0; I--) ASy += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASy);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        ASy = 1.0 / SLACK;
        E *= ScaleFactor;
    }

/* Compute infinity-norm of T for O'Leary's estimate. */
    for (MaxY = 0.0, I = Size; I > 0; I--)
        if (MaxY < ABS(T[I])) MaxY = ABS(T[I]);

/*
 * Part 2.  A* z = y where the * represents the transpose.
 * Recall that A = LU implies A* = U* L*.
 */

/* Forward elimination, U* v = y. */
    for (I = 1; I <= Size; I++)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   T[pElement->Col] -= T[I] * pElement->Real;
            pElement = pElement->NextInRow;
        }
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains v, and scale ||T|| to 1/SLACK. */
    for (ASv = 0.0, I = Size; I > 0; I--) ASv += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASv);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        ASy *= ScaleFactor;
    }

/* Backward Substitution, L* z = v. */
    for (I = Size; I >= 1; I--)
    {   pPivot = Matrix->Diag[I];
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   T[I] -= pElement->Real * T[pElement->Row];
            pElement = pElement->NextInCol;
        }
        T[I] *= pPivot->Real;
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains z. */
    for (ASz = 0.0, I = Size; I > 0; I--) ASz += ABS(T[I]);

#if NOT spCOMPLEX
    FREE( Tm );
#endif

    Linpack = ASy / ASz;
    OLeary = E / MaxY;
    InvNormOfInverse = MIN( Linpack, OLeary );
    return InvNormOfInverse / NormOfMatrix;
#endif /* REAL */
}


#if spCOMPLEX
static double ComplexCondition( Matrix, NormOfMatrix, pError )
MatrixPtr Matrix;
double NormOfMatrix;
int *pError;
{
register ElementPtr pElement;
register ComplexVector T, Tm;
register int I, K, Row;
ElementPtr pPivot;
int Size;
double E, Em, ASp, ASm, ASw, ASy, ASv, ASz, MaxY, ScaleFactor;
double Linpack, OLeary, InvNormOfInverse;
ComplexNumber Wp, Wm;

/* Begin `ComplexCondition'. */

    Size = Matrix->Size;
    T = (ComplexVector)Matrix->Intermediate;
    Tm = ALLOC( ComplexNumber, Size+1 );
    if (Tm == NULL)
    {   *pError = spNO_MEMORY;
        return 0.0;
    }
    for (I = Size; I > 0; I--) T[I].Real = T[I].Imag = 0.0;

/*
 * Part 1.  Ay = e.
 * Solve Ay = LUy = e where e consists of +1 and -1 terms with the sign
 * chosen to maximize the size of w in Lw = e.  Since the terms in w can
 * get very large, scaling is used to avoid overflow.
 */

/* Forward elimination. Solves Lw = e while choosing e. */
    E = 1.0;
    for (I = 1; I <= Size; I++)
    {   pPivot = Matrix->Diag[I];
        if (T[I].Real < 0.0) Em = -E; else Em = E;
        Wm = T[I];
        Wm.Real += Em;
        ASm = CMPLX_1_NORM( Wm );
        CMPLX_MULT_ASSIGN( Wm, *pPivot );
        if (CMPLX_1_NORM(Wm) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(Wm) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            E *= ScaleFactor;
            Em *= ScaleFactor;
            ASm *= ScaleFactor;
            SCLR_MULT_ASSIGN( Wm, ScaleFactor );
        }
        Wp = T[I];
        Wp.Real -= Em;
        ASp = CMPLX_1_NORM( Wp );
        CMPLX_MULT_ASSIGN( Wp, *pPivot );

/* Update T for both values of W, minus value is placed in Tm. */
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   Row = pElement->Row;
            /* Cmplx expr: Tm[Row] = T[Row] - (Wp * *pElement). */
            CMPLX_MULT_SUBT( Tm[Row], Wm, *pElement, T[Row] );
            /* Cmplx expr: T[Row] -= Wp * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[Row], Wm, *pElement );
            ASp += CMPLX_1_NORM(T[Row]);
            ASm += CMPLX_1_NORM(Tm[Row]);
            pElement = pElement->NextInCol;
        }

/* If minus value causes more growth, overwrite T with its values. */
        if (ASm > ASp)
        {   T[I] = Wm;
            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {   T[pElement->Row] = Tm[pElement->Row];
                pElement = pElement->NextInCol;
            }
        }
        else T[I] = Wp;
    }

/* Compute 1-norm of T, which now contains w, and scale ||T|| to 1/SLACK. */
    for (ASw = 0.0, I = Size; I > 0; I--) ASw += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASw);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        E *= ScaleFactor;
    }

/* Backward Substitution. Solves Uy = w.*/
    for (I = Size; I >= 1; I--)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   /* Cmplx expr: T[I] -= T[pElement->Col] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[I], T[pElement->Col], *pElement );
            pElement = pElement->NextInRow;
        }
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            E *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains y, and scale ||T|| to 1/SLACK. */
    for (ASy = 0.0, I = Size; I > 0; I--) ASy += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASy);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        ASy = 1.0 / SLACK;
        E *= ScaleFactor;
    }

/* Compute infinity-norm of T for O'Leary's estimate. */
    for (MaxY = 0.0, I = Size; I > 0; I--)
        if (MaxY < CMPLX_1_NORM(T[I])) MaxY = CMPLX_1_NORM(T[I]);

/*
 * Part 2.  A* z = y where the * represents the transpose.
 * Recall that A = LU implies A* = U* L*.
 */

/* Forward elimination, U* v = y. */
    for (I = 1; I <= Size; I++)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   /* Cmplx expr: T[pElement->Col] -= T[I] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[pElement->Col], T[I], *pElement );
            pElement = pElement->NextInRow;
        }
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains v, and scale ||T|| to 1/SLACK. */
    for (ASv = 0.0, I = Size; I > 0; I--) ASv += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASv);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        ASy *= ScaleFactor;
    }

/* Backward Substitution, L* z = v. */
    for (I = Size; I >= 1; I--)
    {   pPivot = Matrix->Diag[I];
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   /* Cmplx expr: T[I] -= T[pElement->Row] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[I], T[pElement->Row], *pElement );
            pElement = pElement->NextInCol;
        }
        CMPLX_MULT_ASSIGN( T[I], *pPivot );
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains z. */
    for (ASz = 0.0, I = Size; I > 0; I--) ASz += CMPLX_1_NORM(T[I]);

    FREE( Tm );

    Linpack = ASy / ASz;
    OLeary = E / MaxY;
    InvNormOfInverse = MIN( Linpack, OLeary );
    return InvNormOfInverse / NormOfMatrix;
}
#endif /* spCOMPLEX */


double spNorm( eMatrix )
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register int I;
double Max = 0.0, AbsRowSum;

/* Begin `spNorm'. */
    ASSERT( IS_SPARSE(Matrix) AND NOT IS_FACTORED(Matrix) );
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

/* Compute row sums. */
#if REAL
    if (NOT Matrix->Complex)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            AbsRowSum = 0.0;
            while (pElement != NULL)
            {   AbsRowSum += ABS( pElement->Real );
                pElement = pElement->NextInRow;
            }
            if (Max < AbsRowSum) Max = AbsRowSum;
        }
    }
#endif
#if spCOMPLEX
    if (Matrix->Complex)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            AbsRowSum = 0.0;
            while (pElement != NULL)
            {   AbsRowSum += CMPLX_1_NORM( *pElement );
                pElement = pElement->NextInRow;
            }
            if (Max < AbsRowSum) Max = AbsRowSum;
        }
    }
#endif
    return Max;
}
#endif /* CONDITION */


#if STABILITY
double spLargestElement( eMatrix )
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I;
double Mag, AbsColSum, Max = 0.0, MaxRow = 0.0, MaxCol = 0.0;
double Pivot;
ComplexNumber cPivot;
register ElementPtr pElement, pDiag;

/* Begin `spLargestElement'. */
    ASSERT( IS_SPARSE(Matrix) );

#if REAL
    if (Matrix->Factored AND NOT Matrix->Complex)
    {   if (Matrix->Error == spSINGULAR) return 0.0;

/* Find the bound on the size of the largest element over all factorization. */
        for (I = 1; I <= Matrix->Size; I++)
        {   pDiag = Matrix->Diag[I];

/* Lower triangular matrix. */
            Pivot = 1.0 / pDiag->Real;
            Mag = ABS( Pivot );
            if (Mag > MaxRow) MaxRow = Mag;
            pElement = Matrix->FirstInRow[I];
            while (pElement != pDiag)
            {   Mag = ABS( pElement->Real );
                if (Mag > MaxRow) MaxRow = Mag;
                pElement = pElement->NextInRow;
            }

/* Upper triangular matrix. */
            pElement = Matrix->FirstInCol[I];
            AbsColSum = 1.0;  /* Diagonal of U is unity. */
            while (pElement != pDiag)
            {   AbsColSum += ABS( pElement->Real );
                pElement = pElement->NextInCol;
            }
            if (AbsColSum > MaxCol) MaxCol = AbsColSum;
        }
    }
    else if (NOT Matrix->Complex)
    {   for (I = 1; I <= Matrix->Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   Mag = ABS( pElement->Real );
                if (Mag > Max) Max = Mag;
                pElement = pElement->NextInCol;
            }
        }
        return Max;
    }
#endif
#if spCOMPLEX
    if (Matrix->Factored AND Matrix->Complex)
    {   if (Matrix->Error == spSINGULAR) return 0.0;

/* Find the bound on the size of the largest element over all factorization. */
        for (I = 1; I <= Matrix->Size; I++)
        {   pDiag = Matrix->Diag[I];

/* Lower triangular matrix. */
            CMPLX_RECIPROCAL( cPivot, *pDiag );
            Mag = CMPLX_1_NORM( cPivot );
            if (Mag > MaxRow) MaxRow = Mag;
            pElement = Matrix->FirstInRow[I];
            while (pElement != pDiag)
            {   Mag = CMPLX_1_NORM( *pElement );
                if (Mag > MaxRow) MaxRow = Mag;
                pElement = pElement->NextInRow;
            }

/* Upper triangular matrix. */
            pElement = Matrix->FirstInCol[I];
            AbsColSum = 1.0;  /* Diagonal of U is unity. */
            while (pElement != pDiag)
            {   AbsColSum += CMPLX_1_NORM( *pElement );
                pElement = pElement->NextInCol;
            }
            if (AbsColSum > MaxCol) MaxCol = AbsColSum;
        }
    }
    else if (Matrix->Complex)
    {   for (I = 1; I <= Matrix->Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   Mag = CMPLX_1_NORM( *pElement );
                if (Mag > Max) Max = Mag;
                pElement = pElement->NextInCol;
            }
        }
        return Max;
    }
#endif
    return MaxRow * MaxCol;
}


double spRoundoff( eMatrix, Rho )
char *eMatrix;
double Rho;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register int Count, I, MaxCount = 0;
double Reid, Gear;

/* Begin `spRoundoff'. */
    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );

/* Compute Barlow's bound if it is not given. */
    if (Rho < 0.0) Rho = spLargestElement( eMatrix );

/* Find the maximum number of off-diagonals in L if not previously computed. */
    if (Matrix->MaxRowCountInLowerTri < 0)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            Count = 0;
            while (pElement->Col < I)
            {   Count++;
                pElement = pElement->NextInRow;
            }
            if (Count > MaxCount) MaxCount = Count;
        }
        Matrix->MaxRowCountInLowerTri = MaxCount;
    }
    else MaxCount = Matrix->MaxRowCountInLowerTri;

/* Compute error bound. */
    Gear = 1.01*((MaxCount + 1) * Matrix->RelThreshold + 1.0) * SQR(MaxCount);
    Reid = 3.01 * Matrix->Size;

    if (Gear < Reid)
        return (MACHINE_RESOLUTION * Rho * Gear);
    else
        return (MACHINE_RESOLUTION * Rho * Reid);
}
#endif
