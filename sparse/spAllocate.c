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
"@(#)$Header: spAllocate.c,v 1.3 88/06/24 05:00:11 kundert Exp $";
#endif


#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"


char * spCreate (Size, Complex, pError)
    int Size, *pError;
    BOOLEAN Complex;
{
    register unsigned SizePlusOne;
    register MatrixPtr Matrix;
    register int I;
    int AllocatedSize;

    /* Begin `spCreate'. */
    /* Clear error flag. */
    *pError = spOKAY;

    /* Test for valid size. */
    if ((Size < 0) OR (Size == 0 AND NOT EXPANDABLE))
    {
        *pError = spPANIC;
        return NULL;
    }

    /* Test for valid type. */
#if NOT spCOMPLEX
    if (Complex)
    {
        *pError = spPANIC;
        return NULL;
    }
#endif
#if NOT REAL
    if (NOT Complex)
    {
        *pError = spPANIC;
        return NULL;
    }
#endif

    /* Create Matrix. */
    AllocatedSize = MAX (Size, MINIMUM_ALLOCATED_SIZE);
    SizePlusOne = (unsigned) (AllocatedSize + 1);

    if ((Matrix = ALLOC (struct MatrixFrame, 1)) == NULL)
    {
        *pError = spNO_MEMORY;
        return NULL;
    }

    /* Initialize matrix */
    Matrix->ID = SPARSE_ID;
    Matrix->Complex = Complex;
    Matrix->PreviousMatrixWasComplex = Complex;
    Matrix->Factored = NO;
    Matrix->Elements = 0;
    Matrix->Error = *pError;
    Matrix->Fillins = 0;
    Matrix->Reordered = NO;
    Matrix->NeedsOrdering = YES;
    Matrix->NumberOfInterchangesIsOdd = NO;
    Matrix->Partitioned = NO;
    Matrix->RowsLinked = NO;
    Matrix->InternalVectorsAllocated = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->Size = Size;
    Matrix->AllocatedSize = AllocatedSize;
    Matrix->ExtSize = Size;
    Matrix->AllocatedExtSize = AllocatedSize;
    Matrix->CurrentSize = 0;
    Matrix->ExtToIntColMap = NULL;
    Matrix->ExtToIntRowMap = NULL;
    Matrix->IntToExtColMap = NULL;
    Matrix->IntToExtRowMap = NULL;
    Matrix->MarkowitzRow = NULL;
    Matrix->MarkowitzCol = NULL;
    Matrix->MarkowitzProd = NULL;
    Matrix->DoCmplxDirect = NULL;
    Matrix->DoRealDirect = NULL;
    Matrix->Intermediate = NULL;
    Matrix->RelThreshold = DEFAULT_THRESHOLD;
    Matrix->AbsThreshold = 0.0;

    Matrix->TopOfAllocationList = NULL;
    Matrix->RecordsRemaining = 0;
    Matrix->ElementsRemaining = 0;
    Matrix->FillinsRemaining = 0;

    RecordAllocation (Matrix, (char *) Matrix);
    if (Matrix->Error == spNO_MEMORY)
        goto MemoryError;

    /* Take out the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif
    Matrix->TrashCan.Row = 0;
    Matrix->TrashCan.Col = 0;
    Matrix->TrashCan.NextInRow = NULL;
    Matrix->TrashCan.NextInCol = NULL;
#if INITIALIZE
    Matrix->TrashCan.pInitInfo = NULL;
#endif

    /* Allocate space in memory for Diag pointer vector. */
    CALLOC (Matrix->Diag, ElementPtr, SizePlusOne);
    if (Matrix->Diag == NULL)
        goto MemoryError;

    /* Allocate space in memory for FirstInCol pointer vector. */
    CALLOC (Matrix->FirstInCol, ElementPtr, SizePlusOne);
    if (Matrix->FirstInCol == NULL)
        goto MemoryError;

    /* Allocate space in memory for FirstInRow pointer vector. */
    CALLOC (Matrix->FirstInRow, ElementPtr, SizePlusOne);
    if (Matrix->FirstInRow == NULL)
        goto MemoryError;

    /* Allocate space in memory for IntToExtColMap vector. */
    if ((Matrix->IntToExtColMap = ALLOC (int, SizePlusOne)) == NULL)
        goto MemoryError;

    /* Allocate space in memory for IntToExtRowMap vector. */
    if ((Matrix->IntToExtRowMap = ALLOC (int, SizePlusOne)) == NULL)
        goto MemoryError;

    /* Initialize MapIntToExt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {
        Matrix->IntToExtRowMap[I] = I;
        Matrix->IntToExtColMap[I] = I;
    }

#if TRANSLATE
    /* Allocate space in memory for ExtToIntColMap vector. */
    if ((Matrix->ExtToIntColMap = ALLOC (int, SizePlusOne)) == NULL)
        goto MemoryError;

    /* Allocate space in memory for ExtToIntRowMap vector. */
    if ((Matrix->ExtToIntRowMap = ALLOC (int, SizePlusOne)) == NULL)
        goto MemoryError;

    /* Initialize MapExtToInt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {
        Matrix->ExtToIntColMap[I] = -1;
        Matrix->ExtToIntRowMap[I] = -1;
    }
    Matrix->ExtToIntColMap[0] = 0;
    Matrix->ExtToIntRowMap[0] = 0;
#endif

    /* Allocate space for fill-ins and initial set of elements. */
    InitializeElementBlocks (Matrix, SPACE_FOR_ELEMENTS * AllocatedSize,
            SPACE_FOR_FILL_INS * AllocatedSize);
    if (Matrix->Error == spNO_MEMORY)
        goto MemoryError;

    return (char *) Matrix;

MemoryError:

    /* Deallocate matrix and return no pointer to matrix if there is not enough
       memory. */
    *pError = spNO_MEMORY;
    spDestroy ((char *) Matrix);
    return NULL;
}


ElementPtr spcGetElement (Matrix)
    MatrixPtr Matrix;
{
    ElementPtr pElement;

    /* Begin `spcGetElement'. */

    /* Allocate block of MatrixElements if necessary. */
    if (Matrix->ElementsRemaining == 0)
    {
        pElement = ALLOC (struct MatrixElement, ELEMENTS_PER_ALLOCATION);
        RecordAllocation (Matrix, (char *) pElement);
        if (Matrix->Error == spNO_MEMORY)
            return NULL;
        Matrix->ElementsRemaining = ELEMENTS_PER_ALLOCATION;
        Matrix->NextAvailElement = pElement;
    }

    /* Update Element counter and return pointer to Element. */
    Matrix->ElementsRemaining--;
    return Matrix->NextAvailElement++;

}


/* static */ InitializeElementBlocks (Matrix, InitialNumberOfElements,
        NumberOfFillinsExpected)
MatrixPtr Matrix;
int InitialNumberOfElements, NumberOfFillinsExpected;
{
    ElementPtr pElement;

    /* Begin `InitializeElementBlocks'. */

    /* Allocate block of MatrixElements for elements. */
    pElement = ALLOC (struct MatrixElement, InitialNumberOfElements);
    RecordAllocation (Matrix, (char *) pElement);
    if (Matrix->Error == spNO_MEMORY)
        return;
    Matrix->ElementsRemaining = InitialNumberOfElements;
    Matrix->NextAvailElement = pElement;

    /* Allocate block of MatrixElements for fill-ins. */
    pElement = ALLOC (struct MatrixElement, NumberOfFillinsExpected);
    RecordAllocation (Matrix, (char *) pElement);
    if (Matrix->Error == spNO_MEMORY)
        return;
    Matrix->FillinsRemaining = NumberOfFillinsExpected;
    Matrix->NextAvailFillin = pElement;

    /* Allocate a fill-in list structure. */
    Matrix->FirstFillinListNode = ALLOC (struct FillinListNodeStruct, 1);
    RecordAllocation (Matrix, (char *) Matrix->FirstFillinListNode);
    if (Matrix->Error == spNO_MEMORY)
        return;
    Matrix->LastFillinListNode = Matrix->FirstFillinListNode;

    Matrix->FirstFillinListNode->pFillinList = pElement;
    Matrix->FirstFillinListNode->NumberOfFillinsInList = NumberOfFillinsExpected;
    Matrix->FirstFillinListNode->Next = NULL;

    return;
}


ElementPtr spcGetFillin (Matrix)
    MatrixPtr Matrix;
{
    struct FillinListNodeStruct *pListNode;
    ElementPtr pFillins;

    /* Begin `spcGetFillin'. */

#if NOT STRIP OR LINT
    if (Matrix->FillinsRemaining == 0)
        return spcGetElement (Matrix);
#endif
#if STRIP OR LINT

    if (Matrix->FillinsRemaining == 0)
    {
        pListNode = Matrix->LastFillinListNode;

        /* First see if there are any stripped fill-ins left. */
        if (pListNode->Next != NULL)
        {
            Matrix->LastFillinListNode = pListNode = pListNode->Next;
            Matrix->FillinsRemaining = pListNode->NumberOfFillinsInList;
            Matrix->NextAvailFillin = pListNode->pFillinList;
        }
        else
        {
            /* Allocate block of fill-ins. */
            pFillins = ALLOC (struct MatrixElement, ELEMENTS_PER_ALLOCATION);
            RecordAllocation (Matrix, (char *) pFillins);
            if (Matrix->Error == spNO_MEMORY)
                return NULL;
            Matrix->FillinsRemaining = ELEMENTS_PER_ALLOCATION;
            Matrix->NextAvailFillin = pFillins;

            /* Allocate a fill-in list structure. */
            pListNode->Next = ALLOC (struct FillinListNodeStruct, 1);
            RecordAllocation (Matrix, (char *) pListNode->Next);
            if (Matrix->Error == spNO_MEMORY)
                return NULL;
            Matrix->LastFillinListNode = pListNode = pListNode->Next;

            pListNode->pFillinList = pFillins;
            pListNode->NumberOfFillinsInList = ELEMENTS_PER_ALLOCATION;
            pListNode->Next = NULL;
        }
    }
#endif

    /* Update Fill-in counter and return pointer to Fill-in. */
    Matrix->FillinsRemaining--;
    return Matrix->NextAvailFillin++;
}


/* static */ RecordAllocation (Matrix, AllocatedPtr)
    MatrixPtr Matrix;
    char *AllocatedPtr;
{
    /* Begin `RecordAllocation'. */
    /*
     * If Allocated pointer is NULL, assume that malloc returned a NULL pointer,
     * which indicates a spNO_MEMORY error.
     */
    if (AllocatedPtr == NULL)
    {
        Matrix->Error = spNO_MEMORY;
        return;
    }

    /* Allocate block of MatrixElements if necessary. */
    if (Matrix->RecordsRemaining == 0)
    {
        AllocateBlockOfAllocationList (Matrix);
        if (Matrix->Error == spNO_MEMORY)
        {
            FREE (AllocatedPtr);
            return;
        }
    }

    /* Add Allocated pointer to Allocation List. */
    (++Matrix->TopOfAllocationList)->AllocatedPtr = AllocatedPtr;
    Matrix->RecordsRemaining--;
    return;

}


/* static */ AllocateBlockOfAllocationList (Matrix)
    MatrixPtr Matrix;
{
    register int I;
    register AllocationListPtr ListPtr;

    /* Begin `AllocateBlockOfAllocationList'. */
    /* Allocate block of records for allocation list. */
    ListPtr = ALLOC (struct AllocationRecord, (ELEMENTS_PER_ALLOCATION + 1));
    if (ListPtr == NULL)
    {
        Matrix->Error = spNO_MEMORY;
        return;
    }

    /* String entries of allocation list into singly linked list.  List is linked
       such that any record points to the one before it. */

    ListPtr->NextRecord = Matrix->TopOfAllocationList;
    Matrix->TopOfAllocationList = ListPtr;
    ListPtr += ELEMENTS_PER_ALLOCATION;
    for (I = ELEMENTS_PER_ALLOCATION; I > 0; I--)
    {
        ListPtr->NextRecord = ListPtr - 1;
        ListPtr--;
    }

    /* Record allocation of space for allocation list on allocation list. */
    Matrix->TopOfAllocationList->AllocatedPtr = (char *) ListPtr;
    Matrix->RecordsRemaining = ELEMENTS_PER_ALLOCATION;

    return;
}


void spDestroy (eMatrix)
    register char *eMatrix;
{
    MatrixPtr Matrix = (MatrixPtr) eMatrix;
    register AllocationListPtr ListPtr, NextListPtr;

    /* Begin `spDestroy'. */
    ASSERT (IS_SPARSE (Matrix));

    /* Deallocate the vectors that are located in the matrix frame. */
    FREE (Matrix->IntToExtColMap);
    FREE (Matrix->IntToExtRowMap);
    FREE (Matrix->ExtToIntColMap);
    FREE (Matrix->ExtToIntRowMap);
    FREE (Matrix->Diag);
    FREE (Matrix->FirstInRow);
    FREE (Matrix->FirstInCol);
    FREE (Matrix->MarkowitzRow);
    FREE (Matrix->MarkowitzCol);
    FREE (Matrix->MarkowitzProd);
    FREE (Matrix->DoCmplxDirect);
    FREE (Matrix->DoRealDirect);
    FREE (Matrix->Intermediate);

    /* Sequentially step through the list of allocated pointers freeing pointers
     * along the way. */
    ListPtr = Matrix->TopOfAllocationList;
    while (ListPtr != NULL)
    {
        NextListPtr = ListPtr->NextRecord;
        FREE (ListPtr->AllocatedPtr);
        ListPtr = NextListPtr;
    }
    return;
}


int spError (eMatrix)
    char *eMatrix;
{
    /* Begin `spError'. */

    if (eMatrix != NULL)
    {
        ASSERT (((MatrixPtr) eMatrix)->ID == SPARSE_ID);
        return ((MatrixPtr) eMatrix)->Error;
    }
    else
        return spNO_MEMORY;		/* This error may actually be spPANIC,
                                 * no way to tell. */
}


void spWhereSingular (eMatrix, pRow, pCol)
    char *eMatrix;
    int *pRow, *pCol;
{
    MatrixPtr Matrix = (MatrixPtr) eMatrix;

    /* Begin `spWhereSingular'. */
    ASSERT (IS_SPARSE (Matrix));

    if (Matrix->Error == spSINGULAR OR Matrix->Error == spZERO_DIAG)
    {
        *pRow = Matrix->SingularRow;
        *pCol = Matrix->SingularCol;
    }
    else
        *pRow = *pCol = 0;
    return;
}


int spGetSize (eMatrix, External)
    char *eMatrix;
    BOOLEAN External;
{
    MatrixPtr Matrix = (MatrixPtr) eMatrix;

    /* Begin `spGetSize'. */
    ASSERT (IS_SPARSE (Matrix));

#if TRANSLATE
    if (External)
        return Matrix->ExtSize;
    else
        return Matrix->Size;
#else
    return Matrix->Size;
#endif
}


void spSetReal (eMatrix)
    char *eMatrix;
{
    /* Begin `spSetReal'. */

    ASSERT (IS_SPARSE ((MatrixPtr) eMatrix) AND REAL);
    ((MatrixPtr) eMatrix)->Complex = NO;
    return;
}


void spSetComplex (eMatrix)
    char *eMatrix;
{
    /* Begin `spSetComplex'. */

    ASSERT (IS_SPARSE ((MatrixPtr) eMatrix) AND spCOMPLEX);
    ((MatrixPtr) eMatrix)->Complex = YES;
    return;
}


int spFillinCount (eMatrix)
    char *eMatrix;
{
    /* Begin `spFillinCount'. */

    ASSERT (IS_SPARSE ((MatrixPtr) eMatrix));
    return ((MatrixPtr) eMatrix)->Fillins;
}


int spElementCount (eMatrix)
    char *eMatrix;
{
    /* Begin `spElementCount'. */

    ASSERT (IS_SPARSE ((MatrixPtr) eMatrix));
    return ((MatrixPtr) eMatrix)->Elements;
}
