
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
        //var = 50.0;
        cje = 3e-12; 
        cjc = 2e-12; 
    }
};

int main(int argc, char *argv[]) {
    NODE(gnd);
    NODE(n1);
    NODE(n2);
    NODE(n3);
    NODE(n4);
    NODE(n5);
    NODE(n6);
    NODE(n7);
    NODE(n8);
    NODE(n9);

    vpulse vin("vin", n1, gnd, 0.0, 1.0, 1.0e-8, 1.0e-9);
    //vdc vin("vin", n1, gnd, 0.0);
    vdc vcc("vcc", n8, gnd, 12);
    vdc vee("vee", n9, gnd, -12);

    qnl q1("q1", n4, n2, n6);
    qnl q2("q2", n5, n3, n6);

    r rs1("rs1", n1, n2, 1e3);
    r rs2("rs2", n3, gnd, 1e3);
    r rc1("rc1", n4, n8, 10e3);
    r rc2("rc2", n5, n8, 10e3);

    qnl q3("q3", n6, n7, n9);
    qnl q4("q4", n7, n7, n9);

    r rbias("rbias", n7, n8, 20e3);

    op::monitor = 1;
    op();

    plot::output = *argv;
    plot::add("n1", "n5");

    dc::monitor = 100;
    kernel::dc(vin, -0.2, 0.2, 0.0001);
    // kernel::dc(vin, 0.0, 10.0, 0.1);

    // tran::monitor = 100;
    // tran(1.0e-10, 2.0e-8, 0.0, 1.0e-9);
}
