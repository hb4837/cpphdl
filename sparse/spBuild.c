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
"@(#)$Header: spBuild.c,v 1.3 88/06/24 05:00:31 kundert Exp $";
#endif


#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"

static void Translate (MatrixPtr, int *, int *);
static EnlargeMatrix (MatrixPtr, register int);
static ExpandTranslationArrays (MatrixPtr, register int);

void spClear (eMatrix)
char *eMatrix;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement;
  register int I;

  /* Begin `spClear'. */
  ASSERT (IS_SPARSE (Matrix));

  /* Clear matrix. */
#if spCOMPLEX
  if (Matrix->PreviousMatrixWasComplex OR Matrix->Complex)
  {
    for (I = Matrix->Size; I > 0; I--)
    {
      pElement = Matrix->FirstInCol[I];
      while (pElement != NULL)
      {
	pElement->Real = 0.0;
	pElement->Imag = 0.0;
	pElement = pElement->NextInCol;
      }
    }
  }
  else
#endif
  {
    for (I = Matrix->Size; I > 0; I--)
    {
      pElement = Matrix->FirstInCol[I];
      while (pElement != NULL)
      {
	pElement->Real = 0.0;
	pElement = pElement->NextInCol;
      }
    }
  }

  /* Empty the trash. */
  Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
  Matrix->TrashCan.Imag = 0.0;
#endif

  Matrix->Error = spOKAY;
  Matrix->Factored = NO;
  Matrix->SingularCol = 0;
  Matrix->SingularRow = 0;
  Matrix->PreviousMatrixWasComplex = Matrix->Complex;
  return;
}


double * spGetElement (eMatrix, Row, Col)
char *eMatrix;
int Row, Col;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  double *pElement;
  ElementPtr spcFindElementInCol ();
  void Translate ();

  /* Begin `spGetElement'. */
  ASSERT (IS_SPARSE (Matrix) AND Row >= 0 AND Col >= 0);

  if ((Row == 0) OR (Col == 0))
    return &Matrix->TrashCan.Real;

#if NOT TRANSLATE
  ASSERT (Matrix->NeedsOrdering);
#endif

#if TRANSLATE
  Translate (Matrix, &Row, &Col);
  if (Matrix->Error == spNO_MEMORY)
    return NULL;
#endif

#if NOT TRANSLATE
#if NOT EXPANDABLE
  ASSERT (Row <= Matrix->Size AND Col <= Matrix->Size);
#endif

#if EXPANDABLE
  /* Re-size Matrix if necessary. */
  if ((Row > Matrix->Size) OR (Col > Matrix->Size))
    EnlargeMatrix (Matrix, MAX (Row, Col));
  if (Matrix->Error == spNO_MEMORY)
    return NULL;
#endif
#endif

  /*
 * The condition part of the following if statement tests to see if the
 * element resides along the diagonal, if it does then it tests to see
 * if the element has been created yet (Diag pointer not NULL).  The
 * pointer to the element is then assigned to Element after it is cast
 * into a pointer to a RealNumber.  This casting makes the pointer into
 * a pointer to Real.  This statement depends on the fact that Real
 * is the first record in the MatrixElement structure.
 */

  if ((Row != Col) OR ((pElement = (double *) Matrix->Diag[Row]) == NULL))
  {
    /*
 * Element does not exist or does not reside along diagonal.  Search
 * column for element.  As in the if statement above, the pointer to the
 * element which is returned by spcFindElementInCol is cast into a
 * pointer to Real, a RealNumber.
 */
    pElement = (double *) spcFindElementInCol (Matrix,
						 &(Matrix->FirstInCol[Col]),
						   Row, Col, YES);
  }
  return pElement;
}


//
// Added by hba, 23-Sep-98
//
double * spFindElement (eMatrix, Row, Col)
char *eMatrix;
int Row, Col;
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  double *pElement;
  ElementPtr spcFindElementInCol ();
  void Translate ();

  /* Begin `spGetElement'. */
  ASSERT (IS_SPARSE (Matrix) AND Row >= 0 AND Col >= 0);

  if ((Row == 0) OR (Col == 0))
    return &Matrix->TrashCan.Real;

#if NOT TRANSLATE
  ASSERT (Matrix->NeedsOrdering);
#endif

#if TRANSLATE
  Translate (Matrix, &Row, &Col);
  if (Matrix->Error == spNO_MEMORY)
    return NULL;
#endif

#if NOT TRANSLATE
#if NOT EXPANDABLE
  ASSERT (Row <= Matrix->Size AND Col <= Matrix->Size);
#endif

#if EXPANDABLE
  /* Re-size Matrix if necessary. */
  if ((Row > Matrix->Size) OR (Col > Matrix->Size))
    EnlargeMatrix (Matrix, MAX (Row, Col));
  if (Matrix->Error == spNO_MEMORY)
    return NULL;
#endif
#endif

  /*
 * The condition part of the following if statement tests to see if the
 * element resides along the diagonal, if it does then it tests to see
 * if the element has been created yet (Diag pointer not NULL).  The
 * pointer to the element is then assigned to Element after it is cast
 * into a pointer to a RealNumber.  This casting makes the pointer into
 * a pointer to Real.  This statement depends on the fact that Real
 * is the first record in the MatrixElement structure.
 */

  if ((Row != Col) OR ((pElement = (double *) Matrix->Diag[Row]) == NULL))
  {
    /*
 * Element does not exist or does not reside along diagonal.  Search
 * column for element.  As in the if statement above, the pointer to the
 * element which is returned by spcFindElementInCol is cast into a
 * pointer to Real, a RealNumber.
 */
    pElement = (double *) spcFindElementInCol (Matrix,
						 &(Matrix->FirstInCol[Col]),
						   Row, Col, NO);
  }
  return pElement;
}


ElementPtr spcFindElementInCol (Matrix, LastAddr, Row, Col, CreateIfMissing)
MatrixPtr Matrix;
register ElementPtr *LastAddr;
register int Row;
int Col;
BOOLEAN CreateIfMissing;
{
  register ElementPtr pElement;
  ElementPtr spcCreateElement ();

  /* Begin `spcFindElementInCol'. */
  pElement = *LastAddr;

  /* Search for element. */
  while (pElement != NULL)
  {
    if (pElement->Row < Row)
    {
      /* Have not reached element yet. */
      LastAddr = &(pElement->NextInCol);
      pElement = pElement->NextInCol;
    }
    else if (pElement->Row == Row)
    {
      /* Reached element. */
      return pElement;
    }
    else
      break;			/* while loop */
  }

  /* Element does not exist and must be created. */
  if (CreateIfMissing)
    return spcCreateElement (Matrix, Row, Col, LastAddr, NO);
  else
    return NULL;
}


#if TRANSLATE
static void Translate (Matrix, Row, Col)
MatrixPtr Matrix;
int *Row, *Col;
{
  register int IntRow, IntCol, ExtRow, ExtCol;

  /* Begin `Translate'. */
  ExtRow = *Row;
  ExtCol = *Col;

  /* Expand translation arrays if necessary. */
  if ((ExtRow > Matrix->AllocatedExtSize) OR
      (ExtCol > Matrix->AllocatedExtSize))
  {
    ExpandTranslationArrays (Matrix, MAX (ExtRow, ExtCol));
    if (Matrix->Error == spNO_MEMORY)
      return;
  }

  /* Set ExtSize if necessary. */
  if ((ExtRow > Matrix->ExtSize) OR (ExtCol > Matrix->ExtSize))
    Matrix->ExtSize = MAX (ExtRow, ExtCol);

  /* Translate external row or node number to internal row or node number. */
  if ((IntRow = Matrix->ExtToIntRowMap[ExtRow]) == -1)
  {
    Matrix->ExtToIntRowMap[ExtRow] = ++Matrix->CurrentSize;
    Matrix->ExtToIntColMap[ExtRow] = Matrix->CurrentSize;
    IntRow = Matrix->CurrentSize;

#if NOT EXPANDABLE
    ASSERT (IntRow <= Matrix->Size);
#endif

#if EXPANDABLE
    /* Re-size Matrix if necessary. */
    if (IntRow > Matrix->Size)
      EnlargeMatrix (Matrix, IntRow);
    if (Matrix->Error == spNO_MEMORY)
      return;
#endif

    Matrix->IntToExtRowMap[IntRow] = ExtRow;
    Matrix->IntToExtColMap[IntRow] = ExtRow;
  }

  /* Translate external column or node number to internal column or node number.*/
  if ((IntCol = Matrix->ExtToIntColMap[ExtCol]) == -1)
  {
    Matrix->ExtToIntRowMap[ExtCol] = ++Matrix->CurrentSize;
    Matrix->ExtToIntColMap[ExtCol] = Matrix->CurrentSize;
    IntCol = Matrix->CurrentSize;

#if NOT EXPANDABLE
    ASSERT (IntCol <= Matrix->Size);
#endif

#if EXPANDABLE
    /* Re-size Matrix if necessary. */
    if (IntCol > Matrix->Size)
      EnlargeMatrix (Matrix, IntCol);
    if (Matrix->Error == spNO_MEMORY)
      return;
#endif

    Matrix->IntToExtRowMap[IntCol] = ExtCol;
    Matrix->IntToExtColMap[IntCol] = ExtCol;
  }

  *Row = IntRow;
  *Col = IntCol;
  return;
}

#endif


#if QUAD_ELEMENT
int spGetAdmittance (Matrix, Node1, Node2, Template)
char *Matrix;
int Node1, Node2;
struct spTemplate *Template;
{

  /* Begin `spGetAdmittance'. */
  Template->Element1 = spGetElement (Matrix, Node1, Node1);
  Template->Element2 = spGetElement (Matrix, Node2, Node2);
  Template->Element3Negated = spGetElement (Matrix, Node2, Node1);
  Template->Element4Negated = spGetElement (Matrix, Node1, Node2);
  if
    ((Template->Element1 == NULL)
     OR (Template->Element2 == NULL)
     OR (Template->Element3Negated == NULL)
     OR (Template->Element4Negated == NULL)
    )
    return spNO_MEMORY;

  if (Node1 == 0)
    SWAP (double *, Template->Element1, Template->Element2);

  return spOKAY;
}

#endif /* QUAD_ELEMENT */


#if QUAD_ELEMENT
int spGetQuad (Matrix, Row1, Row2, Col1, Col2, Template)
char *Matrix;
int Row1, Row2, Col1, Col2;
struct spTemplate *Template;
{
  /* Begin `spGetQuad'. */
  Template->Element1 = spGetElement (Matrix, Row1, Col1);
  Template->Element2 = spGetElement (Matrix, Row2, Col2);
  Template->Element3Negated = spGetElement (Matrix, Row2, Col1);
  Template->Element4Negated = spGetElement (Matrix, Row1, Col2);
  if
    ((Template->Element1 == NULL)
     OR (Template->Element2 == NULL)
     OR (Template->Element3Negated == NULL)
     OR (Template->Element4Negated == NULL)
    )
    return spNO_MEMORY;

  if (Template->Element1 == &((MatrixPtr) Matrix)->TrashCan.Real)
    SWAP (double *, Template->Element1, Template->Element2);

  return spOKAY;
}
#endif /* QUAD_ELEMENT */


#if QUAD_ELEMENT
int spGetOnes (Matrix, Pos, Neg, Eqn, Template)
char *Matrix;
int Pos, Neg, Eqn;
struct spTemplate *Template;
{
  /* Begin `spGetOnes'. */
  Template->Element4Negated = spGetElement (Matrix, Neg, Eqn);
  Template->Element3Negated = spGetElement (Matrix, Eqn, Neg);
  Template->Element2 = spGetElement (Matrix, Pos, Eqn);
  Template->Element1 = spGetElement (Matrix, Eqn, Pos);
  if
    ((Template->Element1 == NULL)
     OR (Template->Element2 == NULL)
     OR (Template->Element3Negated == NULL)
     OR (Template->Element4Negated == NULL)
    )
    return spNO_MEMORY;

  spADD_REAL_QUAD (*Template, 1.0);
  return spOKAY;
}
#endif /* QUAD_ELEMENT */

ElementPtr spcCreateElement (Matrix, Row, Col, LastAddr, Fillin)
MatrixPtr Matrix;
int Row;
register int Col;
register ElementPtr *LastAddr;
BOOLEAN Fillin;
{
  register ElementPtr pElement, pLastElement;
  ElementPtr pCreatedElement, spcGetElement (), spcGetFillin ();

  /* Begin `spcCreateElement'. */

  if (Matrix->RowsLinked)
  {
    /* Row pointers cannot be ignored. */
    if (Fillin)
    {
      pElement = spcGetFillin (Matrix);
      Matrix->Fillins++;
    }
    else
    {
      pElement = spcGetElement (Matrix);
      Matrix->NeedsOrdering = YES;
    }
    if (pElement == NULL)
      return NULL;

    /* If element is on diagonal, store pointer in Diag. */
    if (Row == Col)
      Matrix->Diag[Row] = pElement;

    /* Initialize Element. */
    pCreatedElement = pElement;
    pElement->Row = Row;
    pElement->Col = Col;
    pElement->Real = 0.0;
#if spCOMPLEX
    pElement->Imag = 0.0;
#endif
#if INITIALIZE
    pElement->pInitInfo = NULL;
#endif

    /* Splice element into column. */
    pElement->NextInCol = *LastAddr;
    *LastAddr = pElement;

    /* Search row for proper element position. */
    pElement = Matrix->FirstInRow[Row];
    pLastElement = NULL;
    while (pElement != NULL)
    {
      /* Search for element row position. */
      if (pElement->Col < Col)
      {
	/* Have not reached desired element. */
	pLastElement = pElement;
	pElement = pElement->NextInRow;
      }
      else
	pElement = NULL;
    }

    /* Splice element into row. */
    pElement = pCreatedElement;
    if (pLastElement == NULL)
    {
      /* Element is first in row. */
      pElement->NextInRow = Matrix->FirstInRow[Row];
      Matrix->FirstInRow[Row] = pElement;
    }
    else
      /* Element is not first in row. */
    {
      pElement->NextInRow = pLastElement->NextInRow;
      pLastElement->NextInRow = pElement;
    }

  }
  else
  {
    /*
 * Matrix has not been factored yet.  Thus get element rather than fill-in.
 * Also, row pointers can be ignored.
 */

    /* Allocate memory for Element. */
    pElement = spcGetElement (Matrix);
    if (pElement == NULL)
      return NULL;

    /* If element is on diagonal, store pointer in Diag. */
    if (Row == Col)
      Matrix->Diag[Row] = pElement;

    /* Initialize Element. */
    pCreatedElement = pElement;
    pElement->Row = Row;
#if DEBUG
    pElement->Col = Col;
#endif
    pElement->Real = 0.0;
#if spCOMPLEX
    pElement->Imag = 0.0;
#endif
#if INITIALIZE
    pElement->pInitInfo = NULL;
#endif

    /* Splice element into column. */
    pElement->NextInCol = *LastAddr;
    *LastAddr = pElement;
  }

  Matrix->Elements++;
  return pCreatedElement;
}


spcLinkRows (Matrix) 
MatrixPtr Matrix;
{
  register ElementPtr pElement, *FirstInRowEntry;
  register ArrayOfElementPtrs FirstInRowArray;
  register int Col;

  /* Begin `spcLinkRows'. */
  FirstInRowArray = Matrix->FirstInRow;
  for (Col = Matrix->Size; Col >= 1; Col--)
  {
    /* Generate row links for the elements in the Col'th column. */
    pElement = Matrix->FirstInCol[Col];

    while (pElement != NULL)
    {
      pElement->Col = Col;
      FirstInRowEntry = &FirstInRowArray[pElement->Row];
      pElement->NextInRow = *FirstInRowEntry;
      *FirstInRowEntry = pElement;
      pElement = pElement->NextInCol;
    }
  }
  Matrix->RowsLinked = YES;
  return;
}



static EnlargeMatrix (Matrix, NewSize)
MatrixPtr Matrix;
register int NewSize;
{
  register int I, OldAllocatedSize = Matrix->AllocatedSize;

  /* Begin `EnlargeMatrix'. */
  Matrix->Size = NewSize;

  if (NewSize <= OldAllocatedSize)
    return;

  /* Expand the matrix frame. */
  NewSize = MAX (NewSize, EXPANSION_FACTOR * OldAllocatedSize);
  Matrix->AllocatedSize = NewSize;

  if ((REALLOC (Matrix->IntToExtColMap, int, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }
  if ((REALLOC (Matrix->IntToExtRowMap, int, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }
  if ((REALLOC (Matrix->Diag, ElementPtr, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }
  if ((REALLOC (Matrix->FirstInCol, ElementPtr, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }
  if ((REALLOC (Matrix->FirstInRow, ElementPtr, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }

  /*
 * Destroy the Markowitz and Intermediate vectors, they will be recreated
 * in spOrderAndFactor().
 */
  FREE (Matrix->MarkowitzRow);
  FREE (Matrix->MarkowitzCol);
  FREE (Matrix->MarkowitzProd);
  FREE (Matrix->DoRealDirect);
  FREE (Matrix->DoCmplxDirect);
  FREE (Matrix->Intermediate);
  Matrix->InternalVectorsAllocated = NO;

  /* Initialize the new portion of the vectors. */
  for (I = OldAllocatedSize + 1; I <= NewSize; I++)
  {
    Matrix->IntToExtColMap[I] = I;
    Matrix->IntToExtRowMap[I] = I;
    Matrix->Diag[I] = NULL;
    Matrix->FirstInRow[I] = NULL;
    Matrix->FirstInCol[I] = NULL;
  }

  return;
}


#if TRANSLATE
static ExpandTranslationArrays (Matrix, NewSize)
MatrixPtr Matrix;
register int NewSize;
{
  register int I, OldAllocatedSize = Matrix->AllocatedExtSize;

  /* Begin `ExpandTranslationArrays'. */
  Matrix->ExtSize = NewSize;

  if (NewSize <= OldAllocatedSize)
    return;

  /* Expand the translation arrays ExtToIntRowMap and ExtToIntColMap. */
  NewSize = MAX (NewSize, EXPANSION_FACTOR * OldAllocatedSize);
  Matrix->AllocatedExtSize = NewSize;

  if ((REALLOC (Matrix->ExtToIntRowMap, int, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }
  if ((REALLOC (Matrix->ExtToIntColMap, int, NewSize + 1)) == NULL)
  {
    Matrix->Error = spNO_MEMORY;
    return;
  }

  /* Initialize the new portion of the vectors. */
  for (I = OldAllocatedSize + 1; I <= NewSize; I++)
  {
    Matrix->ExtToIntRowMap[I] = -1;
    Matrix->ExtToIntColMap[I] = -1;
  }

  return;
}
#endif


#if INITIALIZE
void spInstallInitInfo (pElement, pInitInfo)
double *pElement;
char *pInitInfo;
{
  /* Begin `spInstallInitInfo'. */
  ASSERT (pElement != NULL);

  ((ElementPtr) pElement)->pInitInfo = pInitInfo;
}


char * spGetInitInfo (pElement)
double *pElement;
{
  /* Begin `spGetInitInfo'. */
  ASSERT (pElement != NULL);

  return (char *) ((ElementPtr) pElement)->pInitInfo;
}


int spInitialize (eMatrix, pInit)
char *eMatrix;
int (*pInit) ();
{
  MatrixPtr Matrix = (MatrixPtr) eMatrix;
  register ElementPtr pElement;
  int J, Error, Col;

  /* Begin `spInitialize'. */
  ASSERT (IS_SPARSE (Matrix));

#if spCOMPLEX
  /* Clear imaginary part of matrix if matrix is real but was complex. */
  if (Matrix->PreviousMatrixWasComplex AND NOT Matrix->Complex)
  {
    for (J = Matrix->Size; J > 0; J--)
    {
      pElement = Matrix->FirstInCol[J];
      while (pElement != NULL)
      {
	pElement->Imag = 0.0;
	pElement = pElement->NextInCol;
      }
    }
  }
#endif /* spCOMPLEX */

  /* Initialize the matrix. */
  for (J = Matrix->Size; J > 0; J--)
  {
    pElement = Matrix->FirstInCol[J];
    Col = Matrix->IntToExtColMap[J];
    while (pElement != NULL)
    {
      if (pElement->pInitInfo == NULL)
      {
	pElement->Real = 0.0;
#               if spCOMPLEX
	pElement->Imag = 0.0;
#               endif
      }
      else
      {
	Error = (*pInit) ((double *) pElement, pElement->pInitInfo,
			  Matrix->IntToExtRowMap[pElement->Row], Col);
	if (Error)
	{
	  Matrix->Error = spFATAL;
	  return Error;
	}

      }
      pElement = pElement->NextInCol;
    }
  }

  /* Empty the trash. */
  Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
  Matrix->TrashCan.Imag = 0.0;
#endif

  Matrix->Error = spOKAY;
  Matrix->Factored = NO;
  Matrix->SingularCol = 0;
  Matrix->SingularRow = 0;
  Matrix->PreviousMatrixWasComplex = Matrix->Complex;
  return 0;
}

#endif /* INITIALIZE */
