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
"@(#)$Header: spFactor.c,v 1.3 88/06/24 05:01:12 kundert Exp $";
#endif


#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"


static int FactorComplexMatrix(MatrixPtr);
static CreateInternalVectors(MatrixPtr);
static CountMarkowitz(MatrixPtr, register RealVector, int);
static MarkowitzProducts(MatrixPtr, int);
static ElementPtr SearchForPivot(MatrixPtr, int, int);
static ElementPtr SearchForSingleton(MatrixPtr, int);
static ElementPtr QuicklySearchDiagonal(MatrixPtr, int);
static ElementPtr QuicklySearchDiagonal(MatrixPtr, int);
static ElementPtr SearchDiagonal(MatrixPtr, register int);
static ElementPtr SearchEntireMatrix(MatrixPtr, int);
static double FindLargestInCol(register ElementPtr pElement);
static double FindBiggestInColExclude(MatrixPtr, register ElementPtr pElement, int);
static ExchangeRowsAndCols(MatrixPtr, ElementPtr, int);
static ExchangeColElements(MatrixPtr, int, ElementPtr Element1, int, ElementPtr Element2, int);
static ExchangeRowElements(MatrixPtr, int, ElementPtr Element1, int, ElementPtr Element2, int);
static RealRowColElimination(MatrixPtr, register ElementPtr);
static ComplexRowColElimination(MatrixPtr, register ElementPtr);
static UpdateMarkowitzNumbers(MatrixPtr, ElementPtr);
static ElementPtr CreateFillin(MatrixPtr, int, int);
static int MatrixIsSingular(MatrixPtr, int);
static int ZeroPivot(MatrixPtr, int);
static WriteStatus(MatrixPtr, int);


int spOrderAndFactor (eMatrix, RHS, RelThreshold, AbsThreshold, DiagPivoting)
char *eMatrix;
double RHS[], RelThreshold, AbsThreshold;
BOOLEAN DiagPivoting;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  ElementPtr pPivot;
  int Step, Size, ReorderingRequired;
  ElementPtr SearchForPivot ();
  double LargestInCol, FindLargestInCol ();

  /* Begin `spOrderAndFactor'. */
  ASSERT (IS_VALID (Matrix) AND NOT Matrix->Factored);

  Matrix->Error = spOKAY;
  Size = Matrix->Size;
  if (RelThreshold <= 0.0)
    RelThreshold = Matrix->RelThreshold;
  if (RelThreshold > 1.0)
    RelThreshold = Matrix->RelThreshold;
  Matrix->RelThreshold = RelThreshold;
  if (AbsThreshold < 0.0)
    AbsThreshold = Matrix->AbsThreshold;
  Matrix->AbsThreshold = AbsThreshold;
  ReorderingRequired = NO;

  if (NOT Matrix->NeedsOrdering)
  {
    /* Matrix has been factored before and reordering is not required. */
    for (Step = 1; Step <= Size; Step++)
    {
      pPivot = Matrix->Diag[Step];
      LargestInCol = FindLargestInCol (pPivot->NextInCol);
      if ((LargestInCol * RelThreshold < ELEMENT_MAG (pPivot)))
      {
	if (Matrix->Complex)
	  ComplexRowColElimination (Matrix, pPivot);
	else
	  RealRowColElimination (Matrix, pPivot);
      }
      else
      {
	ReorderingRequired = YES;
	break;			/* for loop */
      }
    }
    if (NOT ReorderingRequired)
      goto Done;
    else
    {
      /*
 * A pivot was not large enough to maintain accuracy,
 * so a partial reordering is required.
 */

#if ANNOTATE >= ON_STRANGE_BEHAVIOR
      printf ("Reordering,  Step = %1d\n", Step);
#endif
    }
  }				/* End of if(NOT Matrix->NeedsOrdering) */
  else
  {
    /*
 * This is the first time the matrix has been factored.  These few statements
 * indicate to the rest of the code that a full reodering is required rather
 * than a partial reordering, which occurs during a failure of a fast
 * factorization.
 */
    Step = 1;
    if (NOT Matrix->RowsLinked)
      spcLinkRows (Matrix);
    if (NOT Matrix->InternalVectorsAllocated)
      CreateInternalVectors (Matrix);
    if (Matrix->Error >= spFATAL)
      return Matrix->Error;
  }

  /* Form initial Markowitz products. */
  CountMarkowitz (Matrix, RHS, Step);
  MarkowitzProducts (Matrix, Step);
  Matrix->MaxRowCountInLowerTri = -1;

  /* Perform reordering and factorization. */
  for (; Step <= Size; Step++)
  {
    pPivot = SearchForPivot (Matrix, Step, DiagPivoting);
    if (pPivot == NULL)
      return MatrixIsSingular (Matrix, Step);
    ExchangeRowsAndCols (Matrix, pPivot, Step);

    if (Matrix->Complex)
      ComplexRowColElimination (Matrix, pPivot);
    else
      RealRowColElimination (Matrix, pPivot);

    if (Matrix->Error >= spFATAL)
      return Matrix->Error;
    UpdateMarkowitzNumbers (Matrix, pPivot);

#if ANNOTATE == FULL
    WriteStatus (Matrix, Step);
#endif
  }

Done:
  Matrix->NeedsOrdering = NO;
  Matrix->Reordered = YES;
  Matrix->Factored = YES;

  return Matrix->Error;
}


int spFactor (eMatrix)
char *eMatrix;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement;
  register ElementPtr pColumn;
  register int Step, Size;
  double Mult;

  /* Begin `spFactor'. */
  ASSERT (IS_VALID (Matrix) AND NOT Matrix->Factored);

  if (Matrix->NeedsOrdering)
  {
    return spOrderAndFactor (eMatrix, (RealVector) NULL,
			     0.0, 0.0, DIAG_PIVOTING_AS_DEFAULT);
  }
  if (NOT Matrix->Partitioned)
    spPartition (eMatrix, spDEFAULT_PARTITION);
#if spCOMPLEX
  if (Matrix->Complex)
    return FactorComplexMatrix (Matrix);
#endif

#if REAL
  Size = Matrix->Size;

  if (Matrix->Diag[1]->Real == 0.0)
    return ZeroPivot (Matrix, 1);
  Matrix->Diag[1]->Real = 1.0 / Matrix->Diag[1]->Real;

  /* Start factorization. */
  for (Step = 2; Step <= Size; Step++)
  {
    if (Matrix->DoRealDirect[Step])
    {				/* Update column using direct addressing scatter-gather. */
      register double *Dest = (double *) Matrix->Intermediate;

      /* Scatter. */
      pElement = Matrix->FirstInCol[Step];
      while (pElement != NULL)
      {
	Dest[pElement->Row] = pElement->Real;
	pElement = pElement->NextInCol;
      }

      /* Update column. */
      pColumn = Matrix->FirstInCol[Step];
      while (pColumn->Row < Step)
      {
	pElement = Matrix->Diag[pColumn->Row];
	pColumn->Real = Dest[pColumn->Row] * pElement->Real;
	while ((pElement = pElement->NextInCol) != NULL)
	  Dest[pElement->Row] -= pColumn->Real * pElement->Real;
	pColumn = pColumn->NextInCol;
      }

      /* Gather. */
      pElement = Matrix->Diag[Step]->NextInCol;
      while (pElement != NULL)
      {
	pElement->Real = Dest[pElement->Row];
	pElement = pElement->NextInCol;
      }

      /* Check for singular matrix. */
      if (Dest[Step] == 0.0)
	return ZeroPivot (Matrix, Step);
      Matrix->Diag[Step]->Real = 1.0 / Dest[Step];
    }
    else
    {				/* Update column using indirect addressing scatter-gather. */
      register double **pDest = (double **) Matrix->Intermediate;

      /* Scatter. */
      pElement = Matrix->FirstInCol[Step];
      while (pElement != NULL)
      {
	pDest[pElement->Row] = &pElement->Real;
	pElement = pElement->NextInCol;
      }

      /* Update column. */
      pColumn = Matrix->FirstInCol[Step];
      while (pColumn->Row < Step)
      {
	pElement = Matrix->Diag[pColumn->Row];
	Mult = (*pDest[pColumn->Row] *= pElement->Real);
	while ((pElement = pElement->NextInCol) != NULL)
	  *pDest[pElement->Row] -= Mult * pElement->Real;
	pColumn = pColumn->NextInCol;
      }

      /* Check for singular matrix. */
      if (Matrix->Diag[Step]->Real == 0.0)
	return ZeroPivot (Matrix, Step);
      Matrix->Diag[Step]->Real = 1.0 / Matrix->Diag[Step]->Real;
    }
  }

  Matrix->Factored = YES;
  return (Matrix->Error = spOKAY);
#endif /* REAL */
}


#if spCOMPLEX
static int FactorComplexMatrix (Matrix)
MatrixPtr Matrix;
{
  register ElementPtr pElement;
  register ElementPtr pColumn;
  register int Step, Size;
  ComplexNumber Mult, Pivot;

  /* Begin `FactorComplexMatrix'. */
  ASSERT (Matrix->Complex);

  Size = Matrix->Size;
  pElement = Matrix->Diag[1];
  if (ELEMENT_MAG (pElement) == 0.0)
    return ZeroPivot (Matrix, 1);
  /* Cmplx expr: *pPivot = 1.0 / *pPivot. */
  CMPLX_RECIPROCAL (*pElement, *pElement);

  /* Start factorization. */
  for (Step = 2; Step <= Size; Step++)
  {
    if (Matrix->DoCmplxDirect[Step])
    {				/* Update column using direct addressing scatter-gather. */
      register ComplexNumber *Dest;
      Dest = (ComplexNumber *) Matrix->Intermediate;

      /* Scatter. */
      pElement = Matrix->FirstInCol[Step];
      while (pElement != NULL)
      {
	Dest[pElement->Row] = *(ComplexNumber *) pElement;
	pElement = pElement->NextInCol;
      }

      /* Update column. */
      pColumn = Matrix->FirstInCol[Step];
      while (pColumn->Row < Step)
      {
	pElement = Matrix->Diag[pColumn->Row];
	/* Cmplx expr: Mult = Dest[pColumn->Row] * (1.0 / *pPivot). */
	CMPLX_MULT (Mult, Dest[pColumn->Row], *pElement);
	CMPLX_ASSIGN (*pColumn, Mult);
	while ((pElement = pElement->NextInCol) != NULL)
	{			/* Cmplx expr: Dest[pElement->Row] -= Mult * pElement */
	  CMPLX_MULT_SUBT_ASSIGN (Dest[pElement->Row], Mult, *pElement);
	}
	pColumn = pColumn->NextInCol;
      }

      /* Gather. */
      pElement = Matrix->Diag[Step]->NextInCol;
      while (pElement != NULL)
      {
	*(ComplexNumber *) pElement = Dest[pElement->Row];
	pElement = pElement->NextInCol;
      }

      /* Check for singular matrix. */
      Pivot = Dest[Step];
      if (CMPLX_1_NORM (Pivot) == 0.0)
	return ZeroPivot (Matrix, Step);
      CMPLX_RECIPROCAL (*Matrix->Diag[Step], Pivot);
    }
    else
    {				/* Update column using direct addressing scatter-gather. */
      register ComplexNumber **pDest;
      pDest = (ComplexNumber **) Matrix->Intermediate;

      /* Scatter. */
      pElement = Matrix->FirstInCol[Step];
      while (pElement != NULL)
      {
	pDest[pElement->Row] = (ComplexNumber *) pElement;
	pElement = pElement->NextInCol;
      }

      /* Update column. */
      pColumn = Matrix->FirstInCol[Step];
      while (pColumn->Row < Step)
      {
	pElement = Matrix->Diag[pColumn->Row];
	/* Cmplx expr: Mult = *pDest[pColumn->Row] * (1.0 / *pPivot). */
	CMPLX_MULT (Mult, *pDest[pColumn->Row], *pElement);
	CMPLX_ASSIGN (*pDest[pColumn->Row], Mult);
	while ((pElement = pElement->NextInCol) != NULL)
	{			/* Cmplx expr: *pDest[pElement->Row] -= Mult * pElement */
	  CMPLX_MULT_SUBT_ASSIGN (*pDest[pElement->Row], Mult, *pElement);
	}
	pColumn = pColumn->NextInCol;
      }

      /* Check for singular matrix. */
      pElement = Matrix->Diag[Step];
      if (ELEMENT_MAG (pElement) == 0.0)
	return ZeroPivot (Matrix, Step);
      CMPLX_RECIPROCAL (*pElement, *pElement);
    }
  }

  Matrix->Factored = YES;
  return (Matrix->Error = spOKAY);
}

#endif /* spCOMPLEX */


void spPartition (eMatrix, Mode)
char *eMatrix;
int Mode;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement, pColumn;
  register int Step, Size;
  register int *Nc, *No, *Nm;
  BOOLEAN *DoRealDirect, *DoCmplxDirect;

  /* Begin `spPartition'. */
  ASSERT (IS_SPARSE (Matrix));
  if (Matrix->Partitioned)
    return;
  Size = Matrix->Size;
  DoRealDirect = Matrix->DoRealDirect;
  DoCmplxDirect = Matrix->DoCmplxDirect;
  Matrix->Partitioned = YES;

  /* If partition is specified by the user, this is easy. */
  if (Mode == spDEFAULT_PARTITION)
    Mode = DEFAULT_PARTITION;
  if (Mode == spDIRECT_PARTITION)
  {
    for (Step = 1; Step <= Size; Step++)
#if REAL
      DoRealDirect[Step] = YES;
#endif
#if spCOMPLEX
    DoCmplxDirect[Step] = YES;
#endif
    return;
  }
  else if (Mode == spINDIRECT_PARTITION)
  {
    for (Step = 1; Step <= Size; Step++)
#if REAL
      DoRealDirect[Step] = NO;
#endif
#if spCOMPLEX
    DoCmplxDirect[Step] = NO;
#endif
    return;
  }
  else
    ASSERT (Mode == spAUTO_PARTITION);

  /* Otherwise, count all operations needed in when factoring matrix. */
  Nc = (int *) Matrix->MarkowitzRow;
  No = (int *) Matrix->MarkowitzCol;
  Nm = (int *) Matrix->MarkowitzProd;

  /* Start mock-factorization. */
  for (Step = 1; Step <= Size; Step++)
  {
    Nc[Step] = No[Step] = Nm[Step] = 0;

    pElement = Matrix->FirstInCol[Step];
    while (pElement != NULL)
    {
      Nc[Step]++;
      pElement = pElement->NextInCol;
    }

    pColumn = Matrix->FirstInCol[Step];
    while (pColumn->Row < Step)
    {
      pElement = Matrix->Diag[pColumn->Row];
      Nm[Step]++;
      while ((pElement = pElement->NextInCol) != NULL)
	No[Step]++;
      pColumn = pColumn->NextInCol;
    }
  }

  for (Step = 1; Step <= Size; Step++)
  {
    /*
 * The following are just estimates based on a count on the number of
 * machine instructions used on each machine to perform the various
 * tasks.  It was assumed that each machine instruction required the
 * same amount of time (I don't believe this is true for the VAX, and
 * have no idea if this is true for the 68000 family).  For optimum
 * performance, these numbers should be tuned to the machine.
 *   Nc is the number of nonzero elements in the column.
 *   Nm is the number of multipliers in the column.
 *   No is the number of operations in the inner loop.
 */

#define generic
#ifdef hp9000s300
#if REAL
    DoRealDirect[Step] = (Nm[Step] + No[Step] > 3 * Nc[Step] - 2 * Nm[Step]);
#endif
#if spCOMPLEX
    /* On the hp350, it is never profitable to use direct for complex. */
    DoCmplxDirect[Step] = NO;
#endif
#undef generic
#endif

#ifdef vax
#if REAL
    DoRealDirect[Step] = (Nm[Step] + No[Step] > 3 * Nc[Step] - 2 * Nm[Step]);
#endif
#if spCOMPLEX
    DoCmplxDirect[Step] = (Nm[Step] + No[Step] > 7 * Nc[Step] - 4 * Nm[Step]);
#endif
#undef generic
#endif

#ifdef generic
#if REAL
    DoRealDirect[Step] = (Nm[Step] + No[Step] > 3 * Nc[Step] - 2 * Nm[Step]);
#endif
#if spCOMPLEX
    DoCmplxDirect[Step] = (Nm[Step] + No[Step] > 7 * Nc[Step] - 4 * Nm[Step]);
#endif
#undef generic
#endif
  }

#if (ANNOTATE == FULL)
  {
    int Ops = 0;
    for (Step = 1; Step <= Size; Step++)
      Ops += No[Step];
    printf ("Operation count for inner loop of factorization = %d.\n", Ops);
  }
#endif
  return;
}


static CreateInternalVectors (Matrix)
MatrixPtr Matrix;
{
  int Size;

  /* Begin `CreateInternalVectors'. */
  /* Create Markowitz arrays. */
  Size = Matrix->Size;

  if (Matrix->MarkowitzRow == NULL)
  {
    if ((Matrix->MarkowitzRow = ALLOC (int, Size + 1)) == NULL)
        Matrix->Error = spNO_MEMORY;
  }
  if (Matrix->MarkowitzCol == NULL)
  {
    if ((Matrix->MarkowitzCol = ALLOC (int, Size + 1)) == NULL)
        Matrix->Error = spNO_MEMORY;
  }
  if (Matrix->MarkowitzProd == NULL)
  {
    if ((Matrix->MarkowitzProd = ALLOC (long, Size + 2)) == NULL)
        Matrix->Error = spNO_MEMORY;
  }

  /* Create DoDirect vectors for use in spFactor(). */
#if REAL
  if (Matrix->DoRealDirect == NULL)
  {
    if ((Matrix->DoRealDirect = ALLOC (BOOLEAN, Size + 1)) == NULL)
      Matrix->Error = spNO_MEMORY;
  }
#endif
#if spCOMPLEX
  if (Matrix->DoCmplxDirect == NULL)
  {
    if ((Matrix->DoCmplxDirect = ALLOC (BOOLEAN, Size + 1)) == NULL)
      Matrix->Error = spNO_MEMORY;
  }
#endif

  /* Create Intermediate vectors for use in MatrixSolve. */
#if spCOMPLEX
  if (Matrix->Intermediate == NULL)
  {
    if ((Matrix->Intermediate = ALLOC (double, 2 * (Size + 1))) == NULL)
      Matrix->Error = spNO_MEMORY;
  }
#else
  if (Matrix->Intermediate == NULL)
  {
    if ((Matrix->Intermediate = ALLOC (double, Size + 1)) == NULL)
      Matrix->Error = spNO_MEMORY;
  }
#endif

  if (Matrix->Error != spNO_MEMORY)
    Matrix->InternalVectorsAllocated = YES;
  return;
}


static CountMarkowitz (Matrix, RHS, Step)
MatrixPtr Matrix;
register RealVector RHS;
int Step;
{
  register int Count, I, Size = Matrix->Size;
  register ElementPtr pElement;
  int ExtRow;

  /* Begin `CountMarkowitz'. */

  /* Correct array pointer for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS OR NOT spCOMPLEX
  if (RHS != NULL)
    --RHS;
#else
  if (RHS != NULL)
  {
    if (Matrix->Complex)
      RHS -= 2;
    else
      --RHS;
  }
#endif
#endif

  /* Generate MarkowitzRow Count for each row. */
  for (I = Step; I <= Size; I++)
  {
    /* Set Count to -1 initially to remove count due to pivot element. */
    Count = -1;
    pElement = Matrix->FirstInRow[I];
    while (pElement != NULL AND pElement->Col < Step)
      pElement = pElement->NextInRow;
    while (pElement != NULL)
    {
      Count++;
      pElement = pElement->NextInRow;
    }

    /* Include nonzero elements in the RHS vector. */
    ExtRow = Matrix->IntToExtRowMap[I];

#if spSEPARATED_COMPLEX_VECTORS OR NOT spCOMPLEX
    if (RHS != NULL)
      if (RHS[ExtRow] != 0.0)
	Count++;
#else
    if (RHS != NULL)
    {
      if (Matrix->Complex)
      {
	if ((RHS[2 * ExtRow] != 0.0) OR (RHS[2 * ExtRow + 1] != 0.0))
	  Count++;
      }
      else if (RHS[I] != 0.0)
	Count++;
    }
#endif
    Matrix->MarkowitzRow[I] = Count;
  }

  /* Generate the MarkowitzCol count for each column. */
  for (I = Step; I <= Size; I++)
  {
    /* Set Count to -1 initially to remove count due to pivot element. */
    Count = -1;
    pElement = Matrix->FirstInCol[I];
    while (pElement != NULL AND pElement->Row < Step)
      pElement = pElement->NextInCol;
    while (pElement != NULL)
    {
      Count++;
      pElement = pElement->NextInCol;
    }
    Matrix->MarkowitzCol[I] = Count;
  }
  return;
}


static MarkowitzProducts (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  register int I, *pMarkowitzRow, *pMarkowitzCol;
  register long Product, *pMarkowitzProduct;
  register int Size = Matrix->Size;
  double fProduct;

  /* Begin `MarkowitzProducts'. */
  Matrix->Singletons = 0;

  pMarkowitzProduct = &(Matrix->MarkowitzProd[Step]);
  pMarkowitzRow = &(Matrix->MarkowitzRow[Step]);
  pMarkowitzCol = &(Matrix->MarkowitzCol[Step]);

  for (I = Step; I <= Size; I++)
  {
    /* If chance of overflow, use real numbers. */
    if ((*pMarkowitzRow > LARGEST_SHORT_INTEGER AND * pMarkowitzCol != 0) OR
	(*pMarkowitzCol > LARGEST_SHORT_INTEGER AND * pMarkowitzRow != 0))
    {
      fProduct = (double) (*pMarkowitzRow++) * (double) (*pMarkowitzCol++);
      if (fProduct >= LARGEST_LONG_INTEGER)
	*pMarkowitzProduct++ = LARGEST_LONG_INTEGER;
      else
	*pMarkowitzProduct++ = fProduct;
    }
    else
    {
      Product = *pMarkowitzRow++ * *pMarkowitzCol++;
      if ((*pMarkowitzProduct++ = Product) == 0)
	Matrix->Singletons++;
    }
  }
  return;
}



static ElementPtr SearchForPivot (Matrix, Step, DiagPivoting)
MatrixPtr Matrix;
int Step, DiagPivoting;
{
  register ElementPtr ChosenPivot;
  ElementPtr SearchForSingleton ();
  ElementPtr QuicklySearchDiagonal ();
  ElementPtr SearchDiagonal ();
  ElementPtr SearchEntireMatrix ();

  /* Begin `SearchForPivot'. */

  /* If singletons exist, look for an acceptable one to use as pivot. */
  if (Matrix->Singletons)
  {
    ChosenPivot = SearchForSingleton (Matrix, Step);
    if (ChosenPivot != NULL)
    {
      Matrix->PivotSelectionMethod = 's';
      return ChosenPivot;
    }
  }

#if DIAGONAL_PIVOTING
  if (DiagPivoting)
  {
    /*
 * Either no singletons exist or they weren't acceptable.  Take quick first
 * pass at searching diagonal.  First search for element on diagonal of
 * remaining submatrix with smallest Markowitz product, then check to see
 * if it okay numerically.  If not, QuicklySearchDiagonal fails.
 */
    ChosenPivot = QuicklySearchDiagonal (Matrix, Step);
    if (ChosenPivot != NULL)
    {
      Matrix->PivotSelectionMethod = 'q';
      return ChosenPivot;
    }

    /*
 * Quick search of diagonal failed, carefully search diagonal and check each
 * pivot candidate numerically before even tentatively accepting it.
 */
    ChosenPivot = SearchDiagonal (Matrix, Step);
    if (ChosenPivot != NULL)
    {
      Matrix->PivotSelectionMethod = 'd';
      return ChosenPivot;
    }
  }
#endif /* DIAGONAL_PIVOTING */

  /* No acceptable pivot found yet, search entire matrix. */
  ChosenPivot = SearchEntireMatrix (Matrix, Step);
  Matrix->PivotSelectionMethod = 'e';

  return ChosenPivot;
}


static ElementPtr SearchForSingleton (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  register ElementPtr ChosenPivot;
  register int I;
  register long *pMarkowitzProduct;
  int Singletons;
  double PivotMag, FindBiggestInColExclude ();

  /* Begin `SearchForSingleton'. */
  /* Initialize pointer that is to scan through MarkowitzProduct vector. */
  pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size + 1]);
  Matrix->MarkowitzProd[Matrix->Size + 1] = Matrix->MarkowitzProd[Step];

  /* Decrement the count of available singletons, on the assumption that an
 * acceptable one will be found. */
  Singletons = Matrix->Singletons--;

  /*
 * Assure that following while loop will always terminate, this is just
 * preventive medicine, if things are working right this should never
 * be needed.
 */
  Matrix->MarkowitzProd[Step - 1] = 0;

  while (Singletons-- > 0)
  {
    /* Singletons exist, find them. */

    /*
 * This is tricky.  Am using a pointer to sequentially step through the
 * MarkowitzProduct array.  Search terminates when singleton (Product = 0)
 * is found.  Note that the conditional in the while statement
 * ( *pMarkowitzProduct ) is true as long as the MarkowitzProduct is not
 * equal to zero.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element, the singleton.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

    while (*pMarkowitzProduct--)
    {				/*
             * N bottles of beer on the wall;
             * N bottles of beer.
             * you take one down and pass it around;
             * N-1 bottles of beer on the wall.
             */
    }
    I = pMarkowitzProduct - Matrix->MarkowitzProd + 1;

    /* Assure that I is valid. */
    if (I < Step)
      break;			/* while (Singletons-- > 0) */
    if (I > Matrix->Size)
      I = Step;

    /* Singleton has been found in either/both row or/and column I. */
    if ((ChosenPivot = Matrix->Diag[I]) != NULL)
    {
      /* Singleton lies on the diagonal. */
      PivotMag = ELEMENT_MAG (ChosenPivot);
      if
	(PivotMag > Matrix->AbsThreshold AND
	 PivotMag > Matrix->RelThreshold *
	 FindBiggestInColExclude (Matrix, ChosenPivot, Step)
	)
	return ChosenPivot;
    }
    else
    {
      /* Singleton does not lie on diagonal, find it. */
      if (Matrix->MarkowitzCol[I] == 0)
      {
	ChosenPivot = Matrix->FirstInCol[I];
	while ((ChosenPivot != NULL) AND (ChosenPivot->Row < Step))
	  ChosenPivot = ChosenPivot->NextInCol;
	PivotMag = ELEMENT_MAG (ChosenPivot);
	if
	  (PivotMag > Matrix->AbsThreshold AND
	   PivotMag > Matrix->RelThreshold *
	   FindBiggestInColExclude (Matrix, ChosenPivot,
				    Step)
	  )
	  return ChosenPivot;
	else
	{
	  if (Matrix->MarkowitzRow[I] == 0)
	  {
	    ChosenPivot = Matrix->FirstInRow[I];
	    while ((ChosenPivot != NULL) AND (ChosenPivot->Col < Step))
	      ChosenPivot = ChosenPivot->NextInRow;
	    PivotMag = ELEMENT_MAG (ChosenPivot);
	    if
	      (PivotMag > Matrix->AbsThreshold AND
	       PivotMag > Matrix->RelThreshold *
	       FindBiggestInColExclude (Matrix,
					ChosenPivot,
					Step)
	      )
	      return ChosenPivot;
	  }
	}
      }
      else
      {
	ChosenPivot = Matrix->FirstInRow[I];
	while ((ChosenPivot != NULL) AND (ChosenPivot->Col < Step))
	  ChosenPivot = ChosenPivot->NextInRow;
	PivotMag = ELEMENT_MAG (ChosenPivot);
	if
	  (PivotMag > Matrix->AbsThreshold AND
	   PivotMag > Matrix->RelThreshold *
	   FindBiggestInColExclude (Matrix, ChosenPivot,
				    Step)
	  )
	  return ChosenPivot;
      }
    }
    /* Singleton not acceptable (too small), try another. */
  }				/* end of while(lSingletons>0) */

  /*
 * All singletons were unacceptable.  Restore Matrix->Singletons count.
 * Initial assumption that an acceptable singleton would be found was wrong.
 */
  Matrix->Singletons++;
  return NULL;
}


#if DIAGONAL_PIVOTING
#if MODIFIED_MARKOWITZ
static ElementPtr QuicklySearchDiagonal (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  register long MinMarkowitzProduct, *pMarkowitzProduct;
  register ElementPtr pDiag, pOtherInRow, pOtherInCol;
  int I, NumberOfTies;
  ElementPtr ChosenPivot, TiedElements[MAX_MARKOWITZ_TIES + 1];
  double Magnitude, LargestInCol, Ratio, MaxRatio;
  double LargestOffDiagonal;
  double FindBiggestInColExclude ();

  /* Begin `QuicklySearchDiagonal'. */
  NumberOfTies = -1;
  MinMarkowitzProduct = LARGEST_LONG_INTEGER;
  pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size + 2]);
  Matrix->MarkowitzProd[Matrix->Size + 1] = Matrix->MarkowitzProd[Step];

  /* Assure that following while loop will always terminate. */
  Matrix->MarkowitzProd[Step - 1] = -1;

  /*
 * This is tricky.  Am using a pointer in the inner while loop to
 * sequentially step through the MarkowitzProduct array.  Search
 * terminates when the Markowitz product of zero placed at location
 * Step-1 is found.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element. The outer for loop is infinite, broken by using
 * break.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

  for (;;)			/* Endless for loop. */
  {
    while (MinMarkowitzProduct < *(--pMarkowitzProduct))
    {				/*
             * N bottles of beer on the wall;
             * N bottles of beer.
             * You take one down and pass it around;
             * N-1 bottles of beer on the wall.
             */
    }

    I = pMarkowitzProduct - Matrix->MarkowitzProd;

    /* Assure that I is valid; if I < Step, terminate search. */
    if (I < Step)
      break;			/* Endless for loop */
    if (I > Matrix->Size)
      I = Step;

    if ((pDiag = Matrix->Diag[I]) == NULL)
      continue;			/* Endless for loop */
    if ((Magnitude = ELEMENT_MAG (pDiag)) <= Matrix->AbsThreshold)
      continue;			/* Endless for loop */

    if (*pMarkowitzProduct == 1)
    {
      /* Case where only one element exists in row and column other than diagonal. */

      /* Find off diagonal elements. */
      pOtherInRow = pDiag->NextInRow;
      pOtherInCol = pDiag->NextInCol;
      if (pOtherInRow == NULL AND pOtherInCol == NULL)
      {
	pOtherInRow = Matrix->FirstInRow[I];
	while (pOtherInRow != NULL)
	{
	  if (pOtherInRow->Col >= Step AND pOtherInRow->Col != I)
	    break;
	  pOtherInRow = pOtherInRow->NextInRow;
	}
	pOtherInCol = Matrix->FirstInCol[I];
	while (pOtherInCol != NULL)
	{
	  if (pOtherInCol->Row >= Step AND pOtherInCol->Row != I)
	    break;
	  pOtherInCol = pOtherInCol->NextInCol;
	}
      }

      /* Accept diagonal as pivot if diagonal is larger than off diagonals and the
 * off diagonals are placed symmetricly. */
      if (pOtherInRow != NULL AND pOtherInCol != NULL)
      {
	if (pOtherInRow->Col == pOtherInCol->Row)
	{
	  LargestOffDiagonal = MAX (ELEMENT_MAG (pOtherInRow),
				    ELEMENT_MAG (pOtherInCol));
	  if (Magnitude >= LargestOffDiagonal)
	  {
	    /* Accept pivot, it is unlikely to contribute excess error. */
	    return pDiag;
	  }
	}
      }
    }

    if (*pMarkowitzProduct < MinMarkowitzProduct)
    {
      /* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
      TiedElements[0] = pDiag;
      MinMarkowitzProduct = *pMarkowitzProduct;
      NumberOfTies = 0;
    }
    else
    {
      /* This case handles Markowitz ties. */
      if (NumberOfTies < MAX_MARKOWITZ_TIES)
      {
	TiedElements[++NumberOfTies] = pDiag;
	if (NumberOfTies >= MinMarkowitzProduct * TIES_MULTIPLIER)
	  break;		/* Endless for loop */
      }
    }
  }				/* End of endless for loop. */

  /* Test to see if any element was chosen as a pivot candidate. */
  if (NumberOfTies < 0)
    return NULL;

  /* Determine which of tied elements is best numerically. */
  ChosenPivot = NULL;
  MaxRatio = 1.0 / Matrix->RelThreshold;

  for (I = 0; I <= NumberOfTies; I++)
  {
    pDiag = TiedElements[I];
    Magnitude = ELEMENT_MAG (pDiag);
    LargestInCol = FindBiggestInColExclude (Matrix, pDiag, Step);
    Ratio = LargestInCol / Magnitude;
    if (Ratio < MaxRatio)
    {
      ChosenPivot = pDiag;
      MaxRatio = Ratio;
    }
  }
  return ChosenPivot;
}

#else /* Not MODIFIED_MARKOWITZ */

static ElementPtr QuicklySearchDiagonal (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  register long MinMarkowitzProduct, *pMarkowitzProduct;
  register ElementPtr pDiag;
  int I;
  ElementPtr ChosenPivot, pOtherInRow, pOtherInCol;
  double Magnitude, LargestInCol, LargestOffDiagonal;
  double FindBiggestInColExclude ();

  /* Begin `QuicklySearchDiagonal'. */
  ChosenPivot = NULL;
  MinMarkowitzProduct = LARGEST_LONG_INTEGER;
  pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size + 2]);
  Matrix->MarkowitzProd[Matrix->Size + 1] = Matrix->MarkowitzProd[Step];

  /* Assure that following while loop will always terminate. */
  Matrix->MarkowitzProd[Step - 1] = -1;

  /*
 * This is tricky.  Am using a pointer in the inner while loop to
 * sequentially step through the MarkowitzProduct array.  Search
 * terminates when the Markowitz product of zero placed at location
 * Step-1 is found.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element. The outer for loop is infinite, broken by using
 * break.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

  for (;;)			/* Endless for loop. */
  {
    while (*(--pMarkowitzProduct) >= MinMarkowitzProduct)
    {				/* Just passing through. */
    }

    I = pMarkowitzProduct - Matrix->MarkowitzProd;

    /* Assure that I is valid; if I < Step, terminate search. */
    if (I < Step)
      break;			/* Endless for loop */
    if (I > Matrix->Size)
      I = Step;

    if ((pDiag = Matrix->Diag[I]) == NULL)
      continue;			/* Endless for loop */
    if ((Magnitude = ELEMENT_MAG (pDiag)) <= Matrix->AbsThreshold)
      continue;			/* Endless for loop */

    if (*pMarkowitzProduct == 1)
    {
      /* Case where only one element exists in row and column other than diagonal. */

      /* Find off-diagonal elements. */
      pOtherInRow = pDiag->NextInRow;
      pOtherInCol = pDiag->NextInCol;
      if (pOtherInRow == NULL AND pOtherInCol == NULL)
      {
	pOtherInRow = Matrix->FirstInRow[I];
	while (pOtherInRow != NULL)
	{
	  if (pOtherInRow->Col >= Step AND pOtherInRow->Col != I)
	    break;
	  pOtherInRow = pOtherInRow->NextInRow;
	}
	pOtherInCol = Matrix->FirstInCol[I];
	while (pOtherInCol != NULL)
	{
	  if (pOtherInCol->Row >= Step AND pOtherInCol->Row != I)
	    break;
	  pOtherInCol = pOtherInCol->NextInCol;
	}
      }

      /* Accept diagonal as pivot if diagonal is larger than off-diagonals and the
 * off-diagonals are placed symmetricly. */
      if (pOtherInRow != NULL AND pOtherInCol != NULL)
      {
	if (pOtherInRow->Col == pOtherInCol->Row)
	{
	  LargestOffDiagonal = MAX (ELEMENT_MAG (pOtherInRow),
				    ELEMENT_MAG (pOtherInCol));
	  if (Magnitude >= LargestOffDiagonal)
	  {
	    /* Accept pivot, it is unlikely to contribute excess error. */
	    return pDiag;
	  }
	}
      }
    }

    MinMarkowitzProduct = *pMarkowitzProduct;
    ChosenPivot = pDiag;
  }				/* End of endless for loop. */

  if (ChosenPivot != NULL)
  {
    LargestInCol = FindBiggestInColExclude (Matrix, ChosenPivot, Step);
    if (ELEMENT_MAG (ChosenPivot) <= Matrix->RelThreshold * LargestInCol)
      ChosenPivot = NULL;
  }
  return ChosenPivot;
}

#endif /* Not MODIFIED_MARKOWITZ */


static ElementPtr SearchDiagonal (Matrix, Step)
MatrixPtr Matrix;
register int Step;
{
  register int J;
  register long MinMarkowitzProduct, *pMarkowitzProduct;
  register int I;
  register ElementPtr pDiag;
  int NumberOfTies, Size = Matrix->Size;
  ElementPtr ChosenPivot;
  double Magnitude, Ratio, RatioOfAccepted, LargestInCol;
  double FindBiggestInColExclude ();

  /* Begin `SearchDiagonal'. */
  ChosenPivot = NULL;
  MinMarkowitzProduct = LARGEST_LONG_INTEGER;
  pMarkowitzProduct = &(Matrix->MarkowitzProd[Size + 2]);
  Matrix->MarkowitzProd[Size + 1] = Matrix->MarkowitzProd[Step];

  /* Start search of diagonal. */
  for (J = Size + 1; J > Step; J--)
  {
    if (*(--pMarkowitzProduct) > MinMarkowitzProduct)
      continue;			/* for loop */
    if (J > Matrix->Size)
      I = Step;
    else
      I = J;
    if ((pDiag = Matrix->Diag[I]) == NULL)
      continue;			/* for loop */
    if ((Magnitude = ELEMENT_MAG (pDiag)) <= Matrix->AbsThreshold)
      continue;			/* for loop */

    /* Test to see if diagonal's magnitude is acceptable. */
    LargestInCol = FindBiggestInColExclude (Matrix, pDiag, Step);
    if (Magnitude <= Matrix->RelThreshold * LargestInCol)
      continue;			/* for loop */

    if (*pMarkowitzProduct < MinMarkowitzProduct)
    {
      /* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
      ChosenPivot = pDiag;
      MinMarkowitzProduct = *pMarkowitzProduct;
      RatioOfAccepted = LargestInCol / Magnitude;
      NumberOfTies = 0;
    }
    else
    {
      /* This case handles Markowitz ties. */
      NumberOfTies++;
      Ratio = LargestInCol / Magnitude;
      if (Ratio < RatioOfAccepted)
      {
	ChosenPivot = pDiag;
	RatioOfAccepted = Ratio;
      }
      if (NumberOfTies >= MinMarkowitzProduct * TIES_MULTIPLIER)
	return ChosenPivot;
    }
  }				/* End of for(Step) */
  return ChosenPivot;
}

#endif /* DIAGONAL_PIVOTING */


static ElementPtr SearchEntireMatrix (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  register int I, Size = Matrix->Size;
  register ElementPtr pElement;
  int NumberOfTies;
  long Product, MinMarkowitzProduct;
  ElementPtr ChosenPivot, pLargestElement;
  double Magnitude, LargestElementMag, Ratio, RatioOfAccepted, LargestInCol;
  double FindLargestInCol ();

  /* Begin `SearchEntireMatrix'. */
  ChosenPivot = NULL;
  LargestElementMag = 0.0;
  MinMarkowitzProduct = LARGEST_LONG_INTEGER;

  /* Start search of matrix on column by column basis. */
  for (I = Step; I <= Size; I++)
  {
    pElement = Matrix->FirstInCol[I];

    while (pElement != NULL AND pElement->Row < Step)
      pElement = pElement->NextInCol;

    if ((LargestInCol = FindLargestInCol (pElement)) == 0.0)
      continue;			/* for loop */

    while (pElement != NULL)
    {
      /* Check to see if element is the largest encountered so far.  If so, record
   its magnitude and address. */
      if ((Magnitude = ELEMENT_MAG (pElement)) > LargestElementMag)
      {
	LargestElementMag = Magnitude;
	pLargestElement = pElement;
      }
      /* Calculate element's MarkowitzProduct. */
      Product = Matrix->MarkowitzRow[pElement->Row] *
	Matrix->MarkowitzCol[pElement->Col];

      /* Test to see if element is acceptable as a pivot candidate. */
      if ((Product <= MinMarkowitzProduct) AND
	  (Magnitude > Matrix->RelThreshold * LargestInCol) AND
	  (Magnitude > Matrix->AbsThreshold))
      {
	/* Test to see if element has lowest MarkowitzProduct yet found, or whether it
   is tied with an element found earlier. */
	if (Product < MinMarkowitzProduct)
	{
	  /* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
	  ChosenPivot = pElement;
	  MinMarkowitzProduct = Product;
	  RatioOfAccepted = LargestInCol / Magnitude;
	  NumberOfTies = 0;
	}
	else
	{
	  /* This case handles Markowitz ties. */
	  NumberOfTies++;
	  Ratio = LargestInCol / Magnitude;
	  if (Ratio < RatioOfAccepted)
	  {
	    ChosenPivot = pElement;
	    RatioOfAccepted = Ratio;
	  }
	  if (NumberOfTies >= MinMarkowitzProduct * TIES_MULTIPLIER)
	    return ChosenPivot;
	}
      }
      pElement = pElement->NextInCol;
    }				/* End of while(pElement != NULL) */
  }				/* End of for(Step) */

  if (ChosenPivot != NULL)
    return ChosenPivot;

  if (LargestElementMag == 0.0)
  {
    Matrix->Error = spSINGULAR;
    return NULL;
  }

  Matrix->Error = spSMALL_PIVOT;
  return pLargestElement;
}


static double FindLargestInCol (pElement)
register ElementPtr pElement;
{
  double Magnitude, Largest = 0.0;

  /* Begin `FindLargestInCol'. */
  /* Search column for largest element beginning at Element. */
  while (pElement != NULL)
  {
    if ((Magnitude = ELEMENT_MAG (pElement)) > Largest)
      Largest = Magnitude;
    pElement = pElement->NextInCol;
  }

  return Largest;
}


static double FindBiggestInColExclude (Matrix, pElement, Step)
MatrixPtr Matrix;
register ElementPtr pElement;
register int Step;
{
  register int Row;
  int Col;
  double Largest, Magnitude;

  /* Begin `FindBiggestInColExclude'. */
  Row = pElement->Row;
  Col = pElement->Col;
  pElement = Matrix->FirstInCol[Col];

  /* Travel down column until reduced submatrix is entered. */
  while ((pElement != NULL) AND (pElement->Row < Step))
    pElement = pElement->NextInCol;

  /* Initialize the variable Largest. */
  if (pElement->Row != Row)
    Largest = ELEMENT_MAG (pElement);
  else
    Largest = 0.0;

  /* Search rest of column for largest element, avoiding excluded element. */
  while ((pElement = pElement->NextInCol) != NULL)
  {
    if ((Magnitude = ELEMENT_MAG (pElement)) > Largest)
    {
      if (pElement->Row != Row)
	Largest = Magnitude;
    }
  }

  return Largest;
}


static ExchangeRowsAndCols (Matrix, pPivot, Step)
MatrixPtr Matrix;
ElementPtr pPivot;
register int Step;
{
  register int Row, Col;
  long OldMarkowitzProd_Step, OldMarkowitzProd_Row, OldMarkowitzProd_Col;
  ElementPtr spcFindElementInCol ();

  /* Begin `ExchangeRowsAndCols'. */
  Row = pPivot->Row;
  Col = pPivot->Col;
  Matrix->PivotsOriginalRow = Row;
  Matrix->PivotsOriginalCol = Col;

  if ((Row == Step) AND (Col == Step))
    return;

  /* Exchange rows and columns. */
  if (Row == Col)
  {
    spcRowExchange (Matrix, Step, Row);
    spcColExchange (Matrix, Step, Col);
    SWAP (long, Matrix->MarkowitzProd[Step], Matrix->MarkowitzProd[Row]);
    SWAP (ElementPtr, Matrix->Diag[Row], Matrix->Diag[Step]);
  }
  else
  {

    /* Initialize variables that hold old Markowitz products. */
    OldMarkowitzProd_Step = Matrix->MarkowitzProd[Step];
    OldMarkowitzProd_Row = Matrix->MarkowitzProd[Row];
    OldMarkowitzProd_Col = Matrix->MarkowitzProd[Col];

    /* Exchange rows. */
    if (Row != Step)
    {
      spcRowExchange (Matrix, Step, Row);
      Matrix->NumberOfInterchangesIsOdd =
	NOT Matrix->NumberOfInterchangesIsOdd;
      Matrix->MarkowitzProd[Row] = Matrix->MarkowitzRow[Row] *
	Matrix->MarkowitzCol[Row];

      /* Update singleton count. */
      if ((Matrix->MarkowitzProd[Row] == 0) != (OldMarkowitzProd_Row == 0))
      {
	if (OldMarkowitzProd_Row == 0)
	  Matrix->Singletons--;
	else
	  Matrix->Singletons++;
      }
    }

    /* Exchange columns. */
    if (Col != Step)
    {
      spcColExchange (Matrix, Step, Col);
      Matrix->NumberOfInterchangesIsOdd =
	NOT Matrix->NumberOfInterchangesIsOdd;
      Matrix->MarkowitzProd[Col] = Matrix->MarkowitzCol[Col] *
	Matrix->MarkowitzRow[Col];

      /* Update singleton count. */
      if ((Matrix->MarkowitzProd[Col] == 0) != (OldMarkowitzProd_Col == 0))
      {
	if (OldMarkowitzProd_Col == 0)
	  Matrix->Singletons--;
	else
	  Matrix->Singletons++;
      }

      Matrix->Diag[Col] = spcFindElementInCol (Matrix,
					       Matrix->FirstInCol + Col,
					       Col, Col, NO);
    }
    if (Row != Step)
    {
      Matrix->Diag[Row] = spcFindElementInCol (Matrix,
					       Matrix->FirstInCol + Row,
					       Row, Row, NO);
    }
    Matrix->Diag[Step] = spcFindElementInCol (Matrix,
					      Matrix->FirstInCol + Step,
					      Step, Step, NO);

    /* Update singleton count. */
    Matrix->MarkowitzProd[Step] = Matrix->MarkowitzCol[Step] *
      Matrix->MarkowitzRow[Step];
    if ((Matrix->MarkowitzProd[Step] == 0) != (OldMarkowitzProd_Step == 0))
    {
      if (OldMarkowitzProd_Step == 0)
	Matrix->Singletons--;
      else
	Matrix->Singletons++;
    }
  }
  return;
}


spcRowExchange (Matrix, Row1, Row2)
MatrixPtr Matrix;
int Row1, Row2;
{
  register ElementPtr Row1Ptr, Row2Ptr;
  int Column;
  ElementPtr Element1, Element2;

  /* Begin `spcRowExchange'. */
  if (Row1 > Row2)
    SWAP (int, Row1, Row2);

  Row1Ptr = Matrix->FirstInRow[Row1];
  Row2Ptr = Matrix->FirstInRow[Row2];
  while (Row1Ptr != NULL OR Row2Ptr != NULL)
  {
    /* Exchange elements in rows while traveling from left to right. */
    if (Row1Ptr == NULL)
    {
      Column = Row2Ptr->Col;
      Element1 = NULL;
      Element2 = Row2Ptr;
      Row2Ptr = Row2Ptr->NextInRow;
    }
    else if (Row2Ptr == NULL)
    {
      Column = Row1Ptr->Col;
      Element1 = Row1Ptr;
      Element2 = NULL;
      Row1Ptr = Row1Ptr->NextInRow;
    }
    else if (Row1Ptr->Col < Row2Ptr->Col)
    {
      Column = Row1Ptr->Col;
      Element1 = Row1Ptr;
      Element2 = NULL;
      Row1Ptr = Row1Ptr->NextInRow;
    }
    else if (Row1Ptr->Col > Row2Ptr->Col)
    {
      Column = Row2Ptr->Col;
      Element1 = NULL;
      Element2 = Row2Ptr;
      Row2Ptr = Row2Ptr->NextInRow;
    }
    else
      /* Row1Ptr->Col == Row2Ptr->Col */
    {
      Column = Row1Ptr->Col;
      Element1 = Row1Ptr;
      Element2 = Row2Ptr;
      Row1Ptr = Row1Ptr->NextInRow;
      Row2Ptr = Row2Ptr->NextInRow;
    }

    ExchangeColElements (Matrix, Row1, Element1, Row2, Element2, Column);
  }				/* end of while(Row1Ptr != NULL OR Row2Ptr != NULL) */

  if (Matrix->InternalVectorsAllocated)
    SWAP (int, Matrix->MarkowitzRow[Row1], Matrix->MarkowitzRow[Row2]);
  SWAP (ElementPtr, Matrix->FirstInRow[Row1], Matrix->FirstInRow[Row2]);
  SWAP (int, Matrix->IntToExtRowMap[Row1], Matrix->IntToExtRowMap[Row2]);
#if TRANSLATE
  Matrix->ExtToIntRowMap[Matrix->IntToExtRowMap[Row1]] = Row1;
  Matrix->ExtToIntRowMap[Matrix->IntToExtRowMap[Row2]] = Row2;
#endif

  return;
}


spcColExchange (Matrix, Col1, Col2)
MatrixPtr Matrix;
int Col1, Col2;
{
  register ElementPtr Col1Ptr, Col2Ptr;
  int Row;
  ElementPtr Element1, Element2;

  /* Begin `spcColExchange'. */
  if (Col1 > Col2)
    SWAP (int, Col1, Col2);

  Col1Ptr = Matrix->FirstInCol[Col1];
  Col2Ptr = Matrix->FirstInCol[Col2];
  while (Col1Ptr != NULL OR Col2Ptr != NULL)
  {
    /* Exchange elements in rows while traveling from top to bottom. */
    if (Col1Ptr == NULL)
    {
      Row = Col2Ptr->Row;
      Element1 = NULL;
      Element2 = Col2Ptr;
      Col2Ptr = Col2Ptr->NextInCol;
    }
    else if (Col2Ptr == NULL)
    {
      Row = Col1Ptr->Row;
      Element1 = Col1Ptr;
      Element2 = NULL;
      Col1Ptr = Col1Ptr->NextInCol;
    }
    else if (Col1Ptr->Row < Col2Ptr->Row)
    {
      Row = Col1Ptr->Row;
      Element1 = Col1Ptr;
      Element2 = NULL;
      Col1Ptr = Col1Ptr->NextInCol;
    }
    else if (Col1Ptr->Row > Col2Ptr->Row)
    {
      Row = Col2Ptr->Row;
      Element1 = NULL;
      Element2 = Col2Ptr;
      Col2Ptr = Col2Ptr->NextInCol;
    }
    else
      /* Col1Ptr->Row == Col2Ptr->Row */
    {
      Row = Col1Ptr->Row;
      Element1 = Col1Ptr;
      Element2 = Col2Ptr;
      Col1Ptr = Col1Ptr->NextInCol;
      Col2Ptr = Col2Ptr->NextInCol;
    }

    ExchangeRowElements (Matrix, Col1, Element1, Col2, Element2, Row);
  }				/* end of while(Col1Ptr != NULL OR Col2Ptr != NULL) */

  if (Matrix->InternalVectorsAllocated)
    SWAP (int, Matrix->MarkowitzCol[Col1], Matrix->MarkowitzCol[Col2]);
  SWAP (ElementPtr, Matrix->FirstInCol[Col1], Matrix->FirstInCol[Col2]);
  SWAP (int, Matrix->IntToExtColMap[Col1], Matrix->IntToExtColMap[Col2]);
#if TRANSLATE
  Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col1]] = Col1;
  Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col2]] = Col2;
#endif

  return;
}


static ExchangeColElements (Matrix, Row1, Element1, Row2, Element2, Column)
MatrixPtr Matrix;
register ElementPtr Element1, Element2;
int Row1, Row2, Column;
{
  ElementPtr *ElementAboveRow1, *ElementAboveRow2;
  ElementPtr ElementBelowRow1, ElementBelowRow2;
  register ElementPtr pElement;

  /* Begin `ExchangeColElements'. */
  /* Search to find the ElementAboveRow1. */
  ElementAboveRow1 = &(Matrix->FirstInCol[Column]);
  pElement = *ElementAboveRow1;
  while (pElement->Row < Row1)
  {
    ElementAboveRow1 = &(pElement->NextInCol);
    pElement = *ElementAboveRow1;
  }
  if (Element1 != NULL)
  {
    ElementBelowRow1 = Element1->NextInCol;
    if (Element2 == NULL)
    {
      /* Element2 does not exist, move Element1 down to Row2. */
      if (ElementBelowRow1 != NULL AND ElementBelowRow1->Row < Row2)
      {
	/* Element1 must be removed from linked list and moved. */
	*ElementAboveRow1 = ElementBelowRow1;

	/* Search column for Row2. */
	pElement = ElementBelowRow1;
	do
	{
	  ElementAboveRow2 = &(pElement->NextInCol);
	  pElement = *ElementAboveRow2;
	}
	while (pElement != NULL AND pElement->Row < Row2);

	/* Place Element1 in Row2. */
	*ElementAboveRow2 = Element1;
	Element1->NextInCol = pElement;
	*ElementAboveRow1 = ElementBelowRow1;
      }
      Element1->Row = Row2;
    }
    else
    {
      /* Element2 does exist, and the two elements must be exchanged. */
      if (ElementBelowRow1->Row == Row2)
      {
	/* Element2 is just below Element1, exchange them. */
	Element1->NextInCol = Element2->NextInCol;
	Element2->NextInCol = Element1;
	*ElementAboveRow1 = Element2;
      }
      else
      {
	/* Element2 is not just below Element1 and must be searched for. */
	pElement = ElementBelowRow1;
	do
	{
	  ElementAboveRow2 = &(pElement->NextInCol);
	  pElement = *ElementAboveRow2;
	}
	while (pElement->Row < Row2);

	ElementBelowRow2 = Element2->NextInCol;

	/* Switch Element1 and Element2. */
	*ElementAboveRow1 = Element2;
	Element2->NextInCol = ElementBelowRow1;
	*ElementAboveRow2 = Element1;
	Element1->NextInCol = ElementBelowRow2;
      }
      Element1->Row = Row2;
      Element2->Row = Row1;
    }
  }
  else
  {
    /* Element1 does not exist. */
    ElementBelowRow1 = pElement;

    /* Find Element2. */
    if (ElementBelowRow1->Row != Row2)
    {
      do
      {
	ElementAboveRow2 = &(pElement->NextInCol);
	pElement = *ElementAboveRow2;
      }
      while (pElement->Row < Row2);

      ElementBelowRow2 = Element2->NextInCol;

      /* Move Element2 to Row1. */
      *ElementAboveRow2 = Element2->NextInCol;
      *ElementAboveRow1 = Element2;
      Element2->NextInCol = ElementBelowRow1;
    }
    Element2->Row = Row1;
  }
  return;
}


static ExchangeRowElements (Matrix, Col1, Element1, Col2, Element2, Row)
MatrixPtr Matrix;
int Col1, Col2, Row;
register ElementPtr Element1, Element2;
{
  ElementPtr *ElementLeftOfCol1, *ElementLeftOfCol2;
  ElementPtr ElementRightOfCol1, ElementRightOfCol2;
  register ElementPtr pElement;

  /* Begin `ExchangeRowElements'. */
  /* Search to find the ElementLeftOfCol1. */
  ElementLeftOfCol1 = &(Matrix->FirstInRow[Row]);
  pElement = *ElementLeftOfCol1;
  while (pElement->Col < Col1)
  {
    ElementLeftOfCol1 = &(pElement->NextInRow);
    pElement = *ElementLeftOfCol1;
  }
  if (Element1 != NULL)
  {
    ElementRightOfCol1 = Element1->NextInRow;
    if (Element2 == NULL)
    {
      /* Element2 does not exist, move Element1 to right to Col2. */
      if (ElementRightOfCol1 != NULL AND ElementRightOfCol1->Col < Col2)
      {
	/* Element1 must be removed from linked list and moved. */
	*ElementLeftOfCol1 = ElementRightOfCol1;

	/* Search Row for Col2. */
	pElement = ElementRightOfCol1;
	do
	{
	  ElementLeftOfCol2 = &(pElement->NextInRow);
	  pElement = *ElementLeftOfCol2;
	}
	while (pElement != NULL AND pElement->Col < Col2);

	/* Place Element1 in Col2. */
	*ElementLeftOfCol2 = Element1;
	Element1->NextInRow = pElement;
	*ElementLeftOfCol1 = ElementRightOfCol1;
      }
      Element1->Col = Col2;
    }
    else
    {
      /* Element2 does exist, and the two elements must be exchanged. */
      if (ElementRightOfCol1->Col == Col2)
      {
	/* Element2 is just right of Element1, exchange them. */
	Element1->NextInRow = Element2->NextInRow;
	Element2->NextInRow = Element1;
	*ElementLeftOfCol1 = Element2;
      }
      else
      {
	/* Element2 is not just right of Element1 and must be searched for. */
	pElement = ElementRightOfCol1;
	do
	{
	  ElementLeftOfCol2 = &(pElement->NextInRow);
	  pElement = *ElementLeftOfCol2;
	}
	while (pElement->Col < Col2);

	ElementRightOfCol2 = Element2->NextInRow;

	/* Switch Element1 and Element2. */
	*ElementLeftOfCol1 = Element2;
	Element2->NextInRow = ElementRightOfCol1;
	*ElementLeftOfCol2 = Element1;
	Element1->NextInRow = ElementRightOfCol2;
      }
      Element1->Col = Col2;
      Element2->Col = Col1;
    }
  }
  else
  {
    /* Element1 does not exist. */
    ElementRightOfCol1 = pElement;

    /* Find Element2. */
    if (ElementRightOfCol1->Col != Col2)
    {
      do
      {
	ElementLeftOfCol2 = &(pElement->NextInRow);
	pElement = *ElementLeftOfCol2;
      }
      while (pElement->Col < Col2);

      ElementRightOfCol2 = Element2->NextInRow;

      /* Move Element2 to Col1. */
      *ElementLeftOfCol2 = Element2->NextInRow;
      *ElementLeftOfCol1 = Element2;
      Element2->NextInRow = ElementRightOfCol1;
    }
    Element2->Col = Col1;
  }
  return;
}


static RealRowColElimination (Matrix, pPivot)
MatrixPtr Matrix;
register ElementPtr pPivot;
{
#if REAL
  register ElementPtr pSub;
  register int Row;
  register ElementPtr pLower, pUpper;
  extern ElementPtr CreateFillin ();

  /* Begin `RealRowColElimination'. */

  /* Test for zero pivot. */
  if (ABS (pPivot->Real) == 0.0)
  {
    (void) MatrixIsSingular (Matrix, pPivot->Row);
    return;
  }
  pPivot->Real = 1.0 / pPivot->Real;

  pUpper = pPivot->NextInRow;
  while (pUpper != NULL)
  {
    /* Calculate upper triangular element. */
    pUpper->Real *= pPivot->Real;

    pSub = pUpper->NextInCol;
    pLower = pPivot->NextInCol;
    while (pLower != NULL)
    {
      Row = pLower->Row;

      /* Find element in row that lines up with current lower triangular element. */
      while (pSub != NULL AND pSub->Row < Row)
	pSub = pSub->NextInCol;

      /* Test to see if desired element was not found, if not, create fill-in. */
      if (pSub == NULL OR pSub->Row > Row)
      {
	pSub = CreateFillin (Matrix, Row, pUpper->Col);
	if (pSub == NULL)
	{
	  Matrix->Error = spNO_MEMORY;
	  return;
	}
      }
      pSub->Real -= pUpper->Real * pLower->Real;
      pSub = pSub->NextInCol;
      pLower = pLower->NextInCol;
    }
    pUpper = pUpper->NextInRow;
  }
  return;
#endif /* REAL */
}


static ComplexRowColElimination (Matrix, pPivot)
MatrixPtr Matrix;
register ElementPtr pPivot;
{
#if spCOMPLEX
  register ElementPtr pSub;
  register int Row;
  register ElementPtr pLower, pUpper;
  ElementPtr CreateFillin ();

  /* Begin `ComplexRowColElimination'. */

  /* Test for zero pivot. */
  if (ELEMENT_MAG (pPivot) == 0.0)
  {
    (void) MatrixIsSingular (Matrix, pPivot->Row);
    return;
  }
  CMPLX_RECIPROCAL (*pPivot, *pPivot);

  pUpper = pPivot->NextInRow;
  while (pUpper != NULL)
  {
    /* Calculate upper triangular element. */
    /* Cmplx expr: *pUpper = *pUpper * (1.0 / *pPivot). */
    CMPLX_MULT_ASSIGN (*pUpper, *pPivot);

    pSub = pUpper->NextInCol;
    pLower = pPivot->NextInCol;
    while (pLower != NULL)
    {
      Row = pLower->Row;

      /* Find element in row that lines up with current lower triangular element. */
      while (pSub != NULL AND pSub->Row < Row)
	pSub = pSub->NextInCol;

      /* Test to see if desired element was not found, if not, create fill-in. */
      if (pSub == NULL OR pSub->Row > Row)
      {
	pSub = CreateFillin (Matrix, Row, pUpper->Col);
	if (pSub == NULL)
	{
	  Matrix->Error = spNO_MEMORY;
	  return;
	}
      }

      /* Cmplx expr: pElement -= *pUpper * pLower. */
      CMPLX_MULT_SUBT_ASSIGN (*pSub, *pUpper, *pLower);
      pSub = pSub->NextInCol;
      pLower = pLower->NextInCol;
    }
    pUpper = pUpper->NextInRow;
  }
  return;
#endif /* spCOMPLEX */
}


static UpdateMarkowitzNumbers (Matrix, pPivot)
MatrixPtr Matrix;
ElementPtr pPivot;
{
  register int Row, Col;
  register ElementPtr ColPtr, RowPtr;
  register int *MarkoRow = Matrix->MarkowitzRow, *MarkoCol = Matrix->MarkowitzCol;
  double Product;

  /* Begin `UpdateMarkowitzNumbers'. */

  /* Update Markowitz numbers. */
  for (ColPtr = pPivot->NextInCol; ColPtr != NULL; ColPtr = ColPtr->NextInCol)
  {
    Row = ColPtr->Row;
    --MarkoRow[Row];

    /* Form Markowitz product while being cautious of overflows. */
    if ((MarkoRow[Row] > LARGEST_SHORT_INTEGER AND MarkoCol[Row] != 0) OR
	(MarkoCol[Row] > LARGEST_SHORT_INTEGER AND MarkoRow[Row] != 0))
    {
      Product = MarkoCol[Row] * MarkoRow[Row];
      if (Product >= LARGEST_LONG_INTEGER)
	Matrix->MarkowitzProd[Row] = LARGEST_LONG_INTEGER;
      else
	Matrix->MarkowitzProd[Row] = Product;
    }
    else
      Matrix->MarkowitzProd[Row] = MarkoRow[Row] * MarkoCol[Row];
    if (MarkoRow[Row] == 0)
      Matrix->Singletons++;
  }

  for (RowPtr = pPivot->NextInRow; RowPtr != NULL; RowPtr = RowPtr->NextInRow)
  {
    Col = RowPtr->Col;
    --MarkoCol[Col];

    /* Form Markowitz product while being cautious of overflows. */
    if ((MarkoRow[Col] > LARGEST_SHORT_INTEGER AND MarkoCol[Col] != 0) OR
	(MarkoCol[Col] > LARGEST_SHORT_INTEGER AND MarkoRow[Col] != 0))
    {
      Product = MarkoCol[Col] * MarkoRow[Col];
      if (Product >= LARGEST_LONG_INTEGER)
	Matrix->MarkowitzProd[Col] = LARGEST_LONG_INTEGER;
      else
	Matrix->MarkowitzProd[Col] = Product;
    }
    else
      Matrix->MarkowitzProd[Col] = MarkoRow[Col] * MarkoCol[Col];
    if ((MarkoCol[Col] == 0) AND (MarkoRow[Col] != 0))
      Matrix->Singletons++;
  }
  return;
}


static ElementPtr CreateFillin (Matrix, Row, Col)
MatrixPtr Matrix;
register int Row;
int Col;
{
  register ElementPtr pElement, *ppElementAbove;
  ElementPtr spcCreateElement ();

  /* Begin `CreateFillin'. */

  /* Find Element above fill-in. */
  ppElementAbove = &Matrix->FirstInCol[Col];
  pElement = *ppElementAbove;
  while (pElement != NULL)
  {
    if (pElement->Row < Row)
    {
      ppElementAbove = &pElement->NextInCol;
      pElement = *ppElementAbove;
    }
    else
      break;			/* while loop */
  }

  /* End of search, create the element. */
  pElement = spcCreateElement (Matrix, Row, Col, ppElementAbove, YES);

  /* Update Markowitz counts and products. */
  Matrix->MarkowitzProd[Row] = ++Matrix->MarkowitzRow[Row] *
    Matrix->MarkowitzCol[Row];
  if ((Matrix->MarkowitzRow[Row] == 1) AND (Matrix->MarkowitzCol[Row] != 0))
    Matrix->Singletons--;
  Matrix->MarkowitzProd[Col] = ++Matrix->MarkowitzCol[Col] *
    Matrix->MarkowitzRow[Col];
  if ((Matrix->MarkowitzRow[Col] != 0) AND (Matrix->MarkowitzCol[Col] == 1))
    Matrix->Singletons--;

  return pElement;
}


static int MatrixIsSingular (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  /* Begin `MatrixIsSingular'. */

  Matrix->SingularRow = Matrix->IntToExtRowMap[Step];
  Matrix->SingularCol = Matrix->IntToExtColMap[Step];
  return (Matrix->Error = spSINGULAR);
}


static int ZeroPivot (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  /* Begin `ZeroPivot'. */

  Matrix->SingularRow = Matrix->IntToExtRowMap[Step];
  Matrix->SingularCol = Matrix->IntToExtColMap[Step];
  return (Matrix->Error = spZERO_DIAG);
}


#if ANNOTATE == FULL
static WriteStatus (Matrix, Step)
MatrixPtr Matrix;
int Step;
{
  int I;

  /* Begin `WriteStatus'. */

  printf ("Step = %1d   ", Step);
  printf ("Pivot found at %1d,%1d using ", Matrix->PivotsOriginalRow,
	  Matrix->PivotsOriginalCol);
  switch (Matrix->PivotSelectionMethod)
  {
  case 's':
    printf ("SearchForSingleton\n");
    break;
  case 'q':
    printf ("QuicklySearchDiagonal\n");
    break;
  case 'd':
    printf ("SearchDiagonal\n");
    break;
  case 'e':
    printf ("SearchEntireMatrix\n");
    break;
  }

  printf ("MarkowitzRow     = ");
  for (I = 1; I <= Matrix->Size; I++)
    printf ("%2d  ", Matrix->MarkowitzRow[I]);
  printf ("\n");

  printf ("MarkowitzCol     = ");
  for (I = 1; I <= Matrix->Size; I++)
    printf ("%2d  ", Matrix->MarkowitzCol[I]);
  printf ("\n");

  printf ("MarkowitzProduct = ");
  for (I = 1; I <= Matrix->Size; I++)
    printf ("%2d  ", Matrix->MarkowitzProd[I]);
  printf ("\n");

  printf ("Singletons = %2d\n", Matrix->Singletons);

  printf ("IntToExtRowMap     = ");
  for (I = 1; I <= Matrix->Size; I++)
    printf ("%2d  ", Matrix->IntToExtRowMap[I]);
  printf ("\n");

  printf ("IntToExtColMap     = ");
  for (I = 1; I <= Matrix->Size; I++)
    printf ("%2d  ", Matrix->IntToExtColMap[I]);
  printf ("\n");

  printf ("ExtToIntRowMap     = ");
  for (I = 1; I <= Matrix->ExtSize; I++)
    printf ("%2d  ", Matrix->ExtToIntRowMap[I]);
  printf ("\n");

  printf ("ExtToIntColMap     = ");
  for (I = 1; I <= Matrix->ExtSize; I++)
    printf ("%2d  ", Matrix->ExtToIntColMap[I]);
  printf ("\n\n");

  /*  spPrint((char *)Matrix, NO, YES);  */

  return;

}

#endif /* ANNOTATE == FULL */
