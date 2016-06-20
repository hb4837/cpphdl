
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

#include "hdl.h"
#include "spice.h"

int main(int argc, char *argv[]) {
    NODE(gnd);
    NODE(n1);

    vdc v1("v1", n1, gnd, 0.0);
    d d1("d1", n1, gnd);
    d1.bv = 5.0;

    plot::output = *argv;
    plot::add("vd(d1)", "id(d1)");
    // dc(v1, -5.4, 1.0, 0.01);
    dc(v1, 0.5, 1.5, 0.0001);
}
