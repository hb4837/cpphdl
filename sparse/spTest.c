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
  e = spGetElement(Matrix,  1,  1); spADD_REAL_ELEMENT(e,  0.0201303);
  e = spGetElement(Matrix,  1,  6); spADD_REAL_ELEMENT(e, -0.0108325);
  e = spGetElement(Matrix,  1, 19); spADD_REAL_ELEMENT(e, -0.00992221);
  e = spGetElement(Matrix,  2,  2); spADD_REAL_ELEMENT(e,  0.00196826);
  e = spGetElement(Matrix,  2,  6); spADD_REAL_ELEMENT(e, -7.00819e-05);
  e = spGetElement(Matrix,  2, 16); spADD_REAL_ELEMENT(e, -0.000149258);
  e = spGetElement(Matrix,  3,  3); spADD_REAL_ELEMENT(e,  0.00233598);
  e = spGetElement(Matrix,  3,  6); spADD_REAL_ELEMENT(e, -6.89538e-05);
  e = spGetElement(Matrix,  3, 13); spADD_REAL_ELEMENT(e, -3.81701e-05);
  e = spGetElement(Matrix,  3, 14); spADD_REAL_ELEMENT(e, -0.000161645);
  e = spGetElement(Matrix,  4,  4); spADD_REAL_ELEMENT(e,  0.00181074);
  e = spGetElement(Matrix,  4,  6); spADD_REAL_ELEMENT(e, -7.11072e-05);
  e = spGetElement(Matrix,  5,  5); spADD_REAL_ELEMENT(e,  0.00193524);
  e = spGetElement(Matrix,  5,  6); spADD_REAL_ELEMENT(e, -8.35603e-05);
  e = spGetElement(Matrix,  5,  8); spADD_REAL_ELEMENT(e, -0.000160605);
  e = spGetElement(Matrix,  6,  1); spADD_REAL_ELEMENT(e, -1.72443e-05);
  e = spGetElement(Matrix,  6,  2); spADD_REAL_ELEMENT(e, -9.41688e-05);
  e = spGetElement(Matrix,  6,  3); spADD_REAL_ELEMENT(e, -4.87603e-05);
  e = spGetElement(Matrix,  6,  4); spADD_REAL_ELEMENT(e, -9.88527e-05);
  e = spGetElement(Matrix,  6,  5); spADD_REAL_ELEMENT(e, -1.01501e-05);
  e = spGetElement(Matrix,  6,  6); spADD_REAL_ELEMENT(e,  0.0237867);
  e = spGetElement(Matrix,  7,  7); spADD_REAL_ELEMENT(e,  0.0020023);
  e = spGetElement(Matrix,  7, 20); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix,  8,  5); spADD_REAL_ELEMENT(e, -6.08416e-05);
  e = spGetElement(Matrix,  8,  8); spADD_REAL_ELEMENT(e,  0.000306666);
  e = spGetElement(Matrix,  8, 11); spADD_REAL_ELEMENT(e,  3.52493e-05);
  e = spGetElement(Matrix,  9,  9); spADD_REAL_ELEMENT(e,  0.00012008);
  e = spGetElement(Matrix, 10, 10); spADD_REAL_ELEMENT(e,  0.00200237);
  e = spGetElement(Matrix, 10, 22); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 11, 11); spADD_REAL_ELEMENT(e,  9.74511e-05);
  e = spGetElement(Matrix, 12, 12); spADD_REAL_ELEMENT(e,  0.000182758);
  e = spGetElement(Matrix, 13,  3); spADD_REAL_ELEMENT(e, -0.000458012);
  e = spGetElement(Matrix, 13, 13); spADD_REAL_ELEMENT(e,  0.00205898);
  e = spGetElement(Matrix, 13, 21); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 14,  3); spADD_REAL_ELEMENT(e, -4.65368e-05);
  e = spGetElement(Matrix, 14, 14); spADD_REAL_ELEMENT(e,  0.000305252);
  e = spGetElement(Matrix, 14, 17); spADD_REAL_ELEMENT(e,  3.7642e-05);
  e = spGetElement(Matrix, 15, 15); spADD_REAL_ELEMENT(e,  0.00012008);
  e = spGetElement(Matrix, 16,  2); spADD_REAL_ELEMENT(e, -0.000156934);
  e = spGetElement(Matrix, 16, 16); spADD_REAL_ELEMENT(e,  0.00217081);
  e = spGetElement(Matrix, 16, 23); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 17, 17); spADD_REAL_ELEMENT(e,  9.5915e-05);
  e = spGetElement(Matrix, 18, 18); spADD_REAL_ELEMENT(e,  0.000182758);
  e = spGetElement(Matrix, 19,  1); spADD_REAL_ELEMENT(e, -0.00413682);
  e = spGetElement(Matrix, 19, 19); spADD_REAL_ELEMENT(e,  0.0189592);
  e = spGetElement(Matrix, 20,  7); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 20, 20); spADD_REAL_ELEMENT(e, 10.0037);
  e = spGetElement(Matrix, 20, 21); spADD_REAL_ELEMENT(e,-10);
  e = spGetElement(Matrix, 21, 13); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 21, 20); spADD_REAL_ELEMENT(e,-10);
  e = spGetElement(Matrix, 21, 21); spADD_REAL_ELEMENT(e, 10.002);
  e = spGetElement(Matrix, 22, 10); spADD_REAL_ELEMENT(e, -0.002);
  e = spGetElement(Matrix, 22, 22); spADD_REAL_ELEMENT(e, 10.0037);
  e = spGetElement(Matrix, 22, 23); spADD_REAL_ELEMENT(e,-10);
  e = spGetElement(Matrix, 23, 16); spADD_REAL_ELEMENT(e,-0.002);
  e = spGetElement(Matrix, 23, 22); spADD_REAL_ELEMENT(e,-10);
  e = spGetElement(Matrix, 23, 23); spADD_REAL_ELEMENT(e, 10.002);
  e = spGetElement(Matrix,  0,  0); spADD_REAL_ELEMENT(e,  0.0);

  spPrint(Matrix, TRUE, FALSE, TRUE);

  RHS[ 1] = -2.57565e-05;
  RHS[ 2] = -4.48769e-07;
  RHS[ 3] = -4.13483e-06;
  RHS[ 4] = -1.01714e-06;
  RHS[ 5] = -1.37342e-05;
  RHS[ 6] = 2.90276e-07;
  RHS[ 7] = 1.51576e-09;
  RHS[ 8] = -1.56322e-07;
  RHS[ 9] = 0;
  RHS[10] = 5.2198e-08;
  RHS[11] = 0;
  RHS[12] = 0;
  RHS[13] = -4.63536e-07;
  RHS[14] = 1.99061e-06;
  RHS[15] = 0;
  RHS[16] = 1.95645e-06;
  RHS[17] = 0;
  RHS[18] = 0;
  RHS[19] = 3.74634e-05;
  RHS[20] = 6.29115e-06;
  RHS[21] = -4.69183e-06;
  RHS[22] = -3.57397e-06;
  RHS[23] = 1.88248e-06;

  spFactor(Matrix);

  spSolve(Matrix, RHS, Solution);

  for (i = 1; i <= 23; i++)
    printf ("%-.9lg\n", (double) Solution[i]);

  exit (0);
}

