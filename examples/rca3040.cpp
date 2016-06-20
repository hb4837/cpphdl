
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
        bf = 80.0;
        vaf = 50.0;
        is = 1.0e-16;
        nf = 1.0;
        nr = 1.0;
        //cje = 3e-12;
        //vje = 0.75;
        //mje = 0.33;
        //cjc = 2e-12;
        //vjc = 0.75;
        //mjc = 0.33;
    }
};

main(int argc, char *argv[]) {
    node gnd;
    node n1("n1");
    node n2("n2");
    node n3("n3");
    node n5("n5");
    node n6("n6");
    node n7("n7");
    node n8("n8");
    node n9("n9");
    node n10("n10");
    node n11("n11");
    node n12("n12");
    node n13("n13");
    node n14("n14");
    node n15("n15");
    node n16("n16");
    node n17("n17");
    node n30("n30");
    node n31("n31");

    vdc vin("vin", n1, gnd, 0.0);
    vdc vcc("vcc", n2, gnd, 15.0);
    vdc vee("vee", n3, gnd, -15.0);

    r rs1("rs1", n30, n1, 1.0e3);
    r rs2("rs2", n31, gnd, 1.0e3);
    r r1("r1", n5, n3, 4.8e3);
    r r2("r2", n6, n3, 4.8e3);
    r r3("r3", n9, n3, 811.0);
    r r4("r4", n8, n3, 2.17e3);
    r r5("r5", n8, gnd, 820.0);
    r r6("r6", n2, n14, 1.32e3);
    r r7("r7", n2, n12, 4.5e3);
    r r8("r8", n2, n15, 1.32e3);
    r r9("r9", n16, gnd, 5.25e3);
    r r10("r10", n17, gnd, 5.25e3);
    qnl q1("q1", n2, n30, n5);
    qnl q2("q2", n2, n31, n6);
    qnl q3("q3", n10, n5, n7);
    qnl q4("q4", n11, n6, n7);
    qnl q5("q5", n14, n12, n10);
    qnl q6("q6", n15, n12, n11);
    qnl q7("q7", n12, n12, n13);
    qnl q8("q8", n13, n13, gnd);
    qnl q9("q9", n7, n8, n9);
    qnl q10("q10", n2, n15, n16);
    qnl q11("q11", n2, n14, n17);

    op::monitor = 1;
    op();

    //plot::output = *argv;
    //plot::add("n16");
    //dc::monitor = 1;
    //dc(vin, -0.2, 0.2, 0.001);
}
