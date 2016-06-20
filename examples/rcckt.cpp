
/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1998,1999
 *  by Bernd Abel
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any usage and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the author is properly credited.  
 *  The author makes no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without expressed
 *  or implied warranty.
 *
 *  The sparse package used in this library is copyrighted by Ken Kundert 
 *  and the University of California, Berkeley. Please refer to the 
 *  copyright notices in the files of the sparse directory.
 *
 *  $Date: 99/01/07$
 *  $Revision: 0.76$
 */

#include "spice.h"

main(int argc, char *argv[]) {
    node gnd;
    node n1("n1");
    node n2("n2");

    vpulse v1("v1", n1, gnd, 0.0, 10.0, 5.0e-9, 1.0e-9);
    r r1("r1", n1, n2, 100.0);
    c c1("c1", n2, gnd, 1.0e-9);

    plot::output = *argv;
    plot::add("n1", "n2", "ipm(c1)");

    tran::order = 2;
    tran::stepsize = variable;
    tran::monitor = 100;
    tran(1.0e-9, 1.0e-6);
}
