
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

class qnl : public q {
public:
    qnl(string s, node& C, node& B, node& E) : q(s, C, B, E, "npn") {
        is = 1e-16; 
        bf = 80.0; 
        br = 1.0; 
        vaf = 50.0;
        cje = 3e-12; 
        cjc = 2e-12; 
    }
};

class qpl : public q {
public:
    qpl(string s, node& C, node& B, node& E) : q(s, C, B, E, "pnp") {
        is = 1e-16; 
        bf = 80.0; 
        vaf = 50.0;
        cje = 3e-12; 
        cjc = 2e-12; 
    }
};

main(int argc, char *argv[]) {
    node gnd;
    node b0("b0", 0.0);
    node b1("b1", 0.0);
    node b2("b2", 0.0);
    node b3("b3", 0.0);
    node c("c", 0.0);

    vdc vce("vce", c, gnd, 0.0);

#ifdef NPN
    i ib0("ib0", b0, gnd, 0.0);
    qnl q0("q0", c, b0, gnd);
    i ib1("ib1", b1, gnd, -1.0e-6);
    qnl q1("q1", c, b1, gnd);
    i ib2("ib2", b2, gnd, -2.0e-6);
    qnl q2("q2", c, b2, gnd);
    i ib3("ib3", b3, gnd, -3.0e-6);
    qnl q3("q3", c, b3, gnd);
#else
    i ib0("ib0", b0, gnd, 0.0);
    qpl q0("q0", c, b0, gnd);
    i ib1("ib1", b1, gnd, 1.0e-6);
    qpl q1("q1", c, b1, gnd);
    i ib2("ib2", b2, gnd, 2.0e-6);
    qpl q2("q2", c, b2, gnd);
    i ib3("ib3", b3, gnd, 3.0e-6);
    qpl q3("q3", c, b3, gnd);
#endif

    plot::output = *argv;
    plot::add("ice(q0)", "ice(q1)", "ice(q2)", "ice(q3)");
    dc(vce, -2.0, 2.0, 0.002);
}
