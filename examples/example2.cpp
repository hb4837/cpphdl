
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

    //i i1("i1", gnd, n1, 1.0);
    //vdc v1("v1", n1, gnd, 2.0);
    vpulse v1("v1", n1, gnd, 0.0, 1.0, 1.0e-7, 1.0e-8, 1.0e-8, 5.0e-7);
    // VCVS
    //e e1("e1", n1, gnd, n2, gnd, 10.0);
    elim elim1("elim1", n1, gnd, n2, gnd, 10.0, -5.0, 5.0);
    // CCVS
    //f f1("f1", n1, gnd, n2, gnd, 0.1);
    // VCCS
    //g g1("g1", n1, gnd, n2, gnd, 0.1);
    // CCVS
    //h h1("h1", n1, gnd, n2, gnd, 0.1);
    //r r2("r2", n2, gnd, 1.0e6);

    plot::output = *argv;
    // plot::add("n1", "n2", "vin(elim1)", "kvin(elim1)");
    plot::add("n1", "n2");

    //op();
    //dc(v1, -1.0, 1.0, 0.001);
    //dc(v1, -1.0, 1.0, 0.00001);
    tran::monitor = 1;
    tran(1.0e-9, 1.0e-6);
}
