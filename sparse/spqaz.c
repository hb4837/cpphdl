#include <stdio.h>
#include <math.h>
#include <ctype.h>
#define spINSIDE_SPARSE
#include "spConfig.h"
#undef spINSIDE_SPARSE
#include "spMatrix.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

main ()
{
  int Error;
  int Last;
  double *e;
  int i;
  static double RHS[100];
  static double Solution[100];
  static char *Matrix;

  Matrix = spCreate (0, 0, &Error);
  if (Matrix == NULL)
  {
    fprintf (stderr, "sparse: insufficient memory available.\n");
    exit(1);
  }

  spClear(Matrix);

  /* mat1, Typical circuit matrix, somewhat ill-conditioned.*/
  e = spGetElement(Matrix,  1,  3); spADD_REAL_ELEMENT(e, 1.0);
  e = spGetElement(Matrix,  1,  4); spADD_REAL_ELEMENT(e, 1.0);

  e = spGetElement(Matrix,  2,  4); spADD_REAL_ELEMENT(e, -1.0);
  e = spGetElement(Matrix,  2,  5); spADD_REAL_ELEMENT(e, 1.0);

  e = spGetElement(Matrix,  3,  1); spADD_REAL_ELEMENT(e, 1.0);

  e = spGetElement(Matrix,  4,  1); spADD_REAL_ELEMENT(e, 0.1);
  e = spGetElement(Matrix,  4,  2); spADD_REAL_ELEMENT(e, -0.1);
  e = spGetElement(Matrix,  4,  4); spADD_REAL_ELEMENT(e, -1.0);

  e = spGetElement(Matrix,  5,  2); spADD_REAL_ELEMENT(e, 0.05);
  e = spGetElement(Matrix,  5,  5); spADD_REAL_ELEMENT(e, -1.0);

  spPrint(Matrix, FALSE, TRUE, TRUE);

  RHS[ 1] = 0.0;
  RHS[ 2] = 0.0;
  RHS[ 3] = -5.0;
  RHS[ 4] = 0.0;
  RHS[ 5] = 0.0;

  spFactor(Matrix);

  spPrint(Matrix, TRUE, TRUE, TRUE);

  spSolve(Matrix, RHS, Solution);

  for (i = 1; i <= 5; i++)
    printf ("%-.9lg\n", (double) Solution[i]);

  exit (0);
}

