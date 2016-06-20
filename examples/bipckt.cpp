
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

class qn : public q {
public:
    qn(string s, node& C, node& B, node& E) : q(s, C, B, E, "npn") {
        bf = 330.0;
        is = 2.0e-13;
    }
};

main(int argc, char *argv[]) {
    node gnd;
    node b("b");
    node e("e");
    node c("c");
    node n4("n4");

    vdc vcc("vcc", n4, gnd, 24.0);
    idc iin("iin", gnd, b, 1.0e-3);
    r r1("r1", b, c, 100.0e3);
    r r2("r2", b, gnd, 6.8e3);
    r r3("r3", c, n4, 1.8e3);
    r r4("r4", e, gnd, 100.0);
    qn q1("q1", c, b, e);

    op::monitor = 1;
    op();
}
