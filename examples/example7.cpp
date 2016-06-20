
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

class qn : public bjt {
public:
    qn(string s, node& C, node& B, node& E) : bjt(s, C, B, E, "npn") {
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
    node n1("n1", 0.0);
    node n2("n2", 0.0);
    node n3("n3", 0.0);
    node n4("n4", 0.0);
    node n5("n5", 0.0);
    node n6("n6", 0.0);

    r r1("r1", n2, gnd, 6.8e3);
    r r2("r2", n3, n2, 100e3);
    r r3("r3", n3, n4, 1.8e3);
    r r4("r5", n5, gnd, 100);
    r ra("ra", n6, gnd, 2.2e3);
    //c c1("c1", n1, n2, 0.47e-6);
    //c c2("c2", n4, n6, 1e-6);
    r rc1("rc1", n1, n2, 1.0e-6);
    r rc2("rc2", n4, n6, 1.0e-6);
    qn q1("q1", n4, n2, n5);
    vdc vh("vh", n3, gnd, 24.0);
    vdc vgen("vgen", n1, gnd, 0.0);
    //vdc vc1("vc1", n2, n1, 1.3832);
    //vdc vc2("vc2", n4, n6, 10.4798);

    //op::monitor = 1;
    //op();

    plot::output = *argv;
    plot::add("n6");
    dc(vgen, -2, 2, 0.01);
}
