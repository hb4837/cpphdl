
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
        is = 2.0e-13; 
        bf = 330; 
        br = 1.0; 
        cje = 0.0; 
        vje = 0.75; 
        mje = 0.33; 
        cjc = 0.0; 
        vjc = 0.75; 
        mjc = 0.33;
    }
};

main(int argc, char *argv[]) {
    node gnd;
    node n1("n1");
    node n2("n2");
    node n3("n3");
    node n4("n4");
    node n5("n5");

    r re("re", n1, n2, 2.2e3);
    vdc vn("vn", n1, gnd, -5.0);
    vdc vp("vp", n5, gnd,  5.0);
    l l1("l1", n4, n5, 5.0e-6);
    c c1("c1", n4, n2, 2.0e-9);
    c c2("c2", n2, n5, 100.0e-12);
    //ipulse i0("i0", n1, n2, 0, 10.0e-6, 0, 0, 0, 25.0e-6);
    qn q1("q1", n4, gnd, n2);

    op::monitor = 1;
    op();

    //plot::output = *argv;
    //plot::add("n4", "n5", "vpm(l1)");

    //tran::monitor = 100;
    //tran(1.0e-9, 10.0e-6);
}
