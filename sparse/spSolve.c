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
"@(#)$Header: spSolve.c,v 1.3 88/06/24 05:02:49 kundert Exp $";
#endif


#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"

static void SolveComplexMatrix(MatrixPtr Matrix, RealVector RHS, RealVector Solution IMAG_VECTORS);
static void SolveComplexTransposedMatrix(MatrixPtr, RealVector RHS, RealVector Solution IMAG_VECTORS);

void spSolve (eMatrix, RHS, Solution IMAG_VECTORS)
char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement;
  register RealVector Intermediate;
  register double Temp;
  register int I, *pExtOrder, Size;
  ElementPtr pPivot;
  void SolveComplexMatrix ();

  /* Begin `spSolve'. */
  ASSERT (IS_VALID (Matrix) AND IS_FACTORED (Matrix));

#if spCOMPLEX
  if (Matrix->Complex)
  {
    SolveComplexMatrix (Matrix, RHS, Solution IMAG_VECTORS);
    return;
  }
#endif

#if REAL
  Intermediate = Matrix->Intermediate;
  Size = Matrix->Size;

  /* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
  --RHS;
  --Solution;
#endif

  /* Initialize Intermediate vector. */
  pExtOrder = &Matrix->IntToExtRowMap[Size];
  for (I = Size; I > 0; I--)
    Intermediate[I] = RHS[*(pExtOrder--)];

  /* Forward elimination. Solves Lc = b.*/
  for (I = 1; I <= Size; I++)
  {
    /* This step of the elimination is skipped if Temp equals zero. */
    if ((Temp = Intermediate[I]) != 0.0)
    {
      pPivot = Matrix->Diag[I];
      Intermediate[I] = (Temp *= pPivot->Real);

      pElement = pPivot->NextInCol;
      while (pElement != NULL)
      {
	Intermediate[pElement->Row] -= Temp * pElement->Real;
	pElement = pElement->NextInCol;
      }
    }
  }

  /* Backward Substitution. Solves Ux = c.*/
  for (I = Size; I > 0; I--)
  {
    Temp = Intermediate[I];
    pElement = Matrix->Diag[I]->NextInRow;
    while (pElement != NULL)
    {
      Temp -= pElement->Real * Intermediate[pElement->Col];
      pElement = pElement->NextInRow;
    }
    Intermediate[I] = Temp;
  }

  /* Unscramble Intermediate vector while placing data in to Solution vector. */
  pExtOrder = &Matrix->IntToExtColMap[Size];
  for (I = Size; I > 0; I--)
    Solution[*(pExtOrder--)] = Intermediate[I];

  return;
#endif /* REAL */
}


#if spCOMPLEX
static void SolveComplexMatrix (Matrix, RHS, Solution IMAG_VECTORS)
MatrixPtr Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
  register ElementPtr pElement;
  register ComplexVector Intermediate;
  register int I, *pExtOrder, Size;
  ElementPtr pPivot;
  register ComplexVector ExtVector;
  ComplexNumber Temp;

  /* Begin `SolveComplexMatrix'. */

  Size = Matrix->Size;
  Intermediate = (ComplexVector) Matrix->Intermediate;

  /* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
  --RHS;
  --iRHS;
  --Solution;
  --iSolution;
#else
  RHS -= 2;
  Solution -= 2;
#endif
#endif

  /* Initialize Intermediate vector. */
  pExtOrder = &Matrix->IntToExtRowMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
  for (I = Size; I > 0; I--)
  {
    Intermediate[I].Real = RHS[*(pExtOrder)];
    Intermediate[I].Imag = iRHS[*(pExtOrder--)];
  }
#else
  ExtVector = (ComplexVector) RHS;
  for (I = Size; I > 0; I--)
    Intermediate[I] = ExtVector[*(pExtOrder--)];
#endif

  /* Forward substitution. Solves Lc = b.*/
  for (I = 1; I <= Size; I++)
  {
    Temp = Intermediate[I];

    /* This step of the substitution is skipped if Temp equals zero. */
    if ((Temp.Real != 0.0) OR (Temp.Imag != 0.0))
    {
      pPivot = Matrix->Diag[I];
      /* Cmplx expr: Temp *= (1.0 / Pivot). */
      CMPLX_MULT_ASSIGN (Temp, *pPivot);
      Intermediate[I] = Temp;
      pElement = pPivot->NextInCol;
      while (pElement != NULL)
      {
	/* Cmplx expr: Intermediate[Element->Row] -= Temp * *Element. */
	CMPLX_MULT_SUBT_ASSIGN (Intermediate[pElement->Row],
				Temp, *pElement);
	pElement = pElement->NextInCol;
      }
    }
  }

  /* Backward Substitution. Solves Ux = c.*/
  for (I = Size; I > 0; I--)
  {
    Temp = Intermediate[I];
    pElement = Matrix->Diag[I]->NextInRow;

    while (pElement != NULL)
    {
      /* Cmplx expr: Temp -= *Element * Intermediate[Element->Col]. */
      CMPLX_MULT_SUBT_ASSIGN (Temp, *pElement, Intermediate[pElement->Col]);
      pElement = pElement->NextInRow;
    }
    Intermediate[I] = Temp;
  }

  /* Unscramble Intermediate vector while placing data in to Solution vector. */
  pExtOrder = &Matrix->IntToExtColMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
  for (I = Size; I > 0; I--)
  {
    Solution[*(pExtOrder)] = Intermediate[I].Real;
    iSolution[*(pExtOrder--)] = Intermediate[I].Imag;
  }
#else
  ExtVector = (ComplexVector) Solution;
  for (I = Size; I > 0; I--)
    ExtVector[*(pExtOrder--)] = Intermediate[I];
#endif

  return;
}

#endif /* spCOMPLEX */


#if TRANSPOSE
void spSolveTransposed (eMatrix, RHS, Solution IMAG_VECTORS)
char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement;
  register RealVector Intermediate;
  register int I, *pExtOrder, Size;
  ElementPtr pPivot;
  double Temp;
  void SolveComplexTransposedMatrix ();

  /* Begin `spSolveTransposed'. */
  ASSERT (IS_VALID (Matrix) AND IS_FACTORED (Matrix));

#if spCOMPLEX
  if (Matrix->Complex)
  {
    SolveComplexTransposedMatrix (Matrix, RHS, Solution IMAG_VECTORS);
    return;
  }
#endif

#if REAL
  Size = Matrix->Size;
  Intermediate = Matrix->Intermediate;

  /* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
  --RHS;
  --Solution;
#endif

  /* Initialize Intermediate vector. */
  pExtOrder = &Matrix->IntToExtColMap[Size];
  for (I = Size; I > 0; I--)
    Intermediate[I] = RHS[*(pExtOrder--)];

  /* Forward elimination. */
  for (I = 1; I <= Size; I++)
  {
    /* This step of the elimination is skipped if Temp equals zero. */
    if ((Temp = Intermediate[I]) != 0.0)
    {
      pElement = Matrix->Diag[I]->NextInRow;
      while (pElement != NULL)
      {
	Intermediate[pElement->Col] -= Temp * pElement->Real;
	pElement = pElement->NextInRow;
      }

    }
  }

  /* Backward Substitution. */
  for (I = Size; I > 0; I--)
  {
    pPivot = Matrix->Diag[I];
    Temp = Intermediate[I];
    pElement = pPivot->NextInCol;
    while (pElement != NULL)
    {
      Temp -= pElement->Real * Intermediate[pElement->Row];
      pElement = pElement->NextInCol;
    }
    Intermediate[I] = Temp * pPivot->Real;
  }

  /* Unscramble Intermediate vector while placing data in to Solution vector. */
  pExtOrder = &Matrix->IntToExtRowMap[Size];
  for (I = Size; I > 0; I--)
    Solution[*(pExtOrder--)] = Intermediate[I];

  return;
#endif /* REAL */
}

#endif /* TRANSPOSE */


#if TRANSPOSE AND spCOMPLEX
static void SolveComplexTransposedMatrix (Matrix, RHS, Solution IMAG_VECTORS)
MatrixPtr Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
  register ElementPtr pElement;
  register ComplexVector Intermediate;
  register int I, *pExtOrder, Size;
  register ComplexVector ExtVector;
  ElementPtr pPivot;
  ComplexNumber Temp;

  /* Begin `SolveComplexTransposedMatrix'. */

  Size = Matrix->Size;
  Intermediate = (ComplexVector) Matrix->Intermediate;

  /* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
  --RHS;
  --iRHS;
  --Solution;
  --iSolution;
#else
  RHS -= 2;
  Solution -= 2;
#endif
#endif

  /* Initialize Intermediate vector. */
  pExtOrder = &Matrix->IntToExtColMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
  for (I = Size; I > 0; I--)
  {
    Intermediate[I].Real = RHS[*(pExtOrder)];
    Intermediate[I].Imag = iRHS[*(pExtOrder--)];
  }
#else
  ExtVector = (ComplexVector) RHS;
  for (I = Size; I > 0; I--)
    Intermediate[I] = ExtVector[*(pExtOrder--)];
#endif

  /* Forward elimination. */
  for (I = 1; I <= Size; I++)
  {
    Temp = Intermediate[I];

    /* This step of the elimination is skipped if Temp equals zero. */
    if ((Temp.Real != 0.0) OR (Temp.Imag != 0.0))
    {
      pElement = Matrix->Diag[I]->NextInRow;
      while (pElement != NULL)
      {
	/* Cmplx expr: Intermediate[Element->Col] -= Temp * *Element. */
	CMPLX_MULT_SUBT_ASSIGN (Intermediate[pElement->Col],
				Temp, *pElement);
	pElement = pElement->NextInRow;
      }
    }
  }

  /* Backward Substitution. */
  for (I = Size; I > 0; I--)
  {
    pPivot = Matrix->Diag[I];
    Temp = Intermediate[I];
    pElement = pPivot->NextInCol;

    while (pElement != NULL)
    {
      /* Cmplx expr: Temp -= Intermediate[Element->Row] * *Element. */
      CMPLX_MULT_SUBT_ASSIGN (Temp, Intermediate[pElement->Row], *pElement);

      pElement = pElement->NextInCol;
    }
    /* Cmplx expr: Intermediate = Temp * (1.0 / *pPivot). */
    CMPLX_MULT (Intermediate[I], Temp, *pPivot);
  }

  /* Unscramble Intermediate vector while placing data in to Solution vector. */
  pExtOrder = &Matrix->IntToExtRowMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
  for (I = Size; I > 0; I--)
  {
    Solution[*(pExtOrder)] = Intermediate[I].Real;
    iSolution[*(pExtOrder--)] = Intermediate[I].Imag;
  }
#else
  ExtVector = (ComplexVector) Solution;
  for (I = Size; I > 0; I--)
    ExtVector[*(pExtOrder--)] = Intermediate[I];
#endif

  return;
}

#endif /* TRANSPOSE AND spCOMPLEX */
