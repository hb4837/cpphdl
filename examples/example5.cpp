
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
    node n1("n1", 0.0);
    node n2("n2", 0.0);

    vpulse v1("v1", n1, gnd, 0.0, 1.0, 0.0, 1.0);
    vco vco1("vco1", n1, gnd, n2, gnd, 5.0, 0.0, 1.0e2);

    plot::output = *argv;
    plot::add("n2", "integ(vco1)");

    op::monitor = 0;
    op();

    tran::monitor = 0;
    tran(1.0e-3, 1.0);
}
