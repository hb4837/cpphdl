
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
        cje = 3e-12;
        cjc = 2e-12;
    }
};

class qpl : public q {
public:
    qpl(string s, node& C, node& B, node& E) : q(s, C, B, E, "pnp") {
        bf = 10.0;
        vaf = 50.0;
        cje = 6e-12;
        cjc = 4e-12;
    }
};

main(int argc, char *argv[]) {
    node gnd;
    node n1("n1");
    node n2("n2");
    node n3("n3");
    node n4("n4");
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
    node n18("n18");
    node n19("n19");
    node n20("n20");
    node n21("n21");
    node n22("n22");
    node n23("n23");
    node n30("n30");
    node n31("n31");

    vdc vin("vin", n1, gnd, 0.0);
    vdc vcc("vcc", n19, gnd, 15.0);
    vdc vee("vee", n20, gnd, -15.0);

    r rs1("rs1", n30, n1, 1.0e3);
    r rs2("rs2", n31, gnd, 1.0e3);
    r rf("rf", n30, n18, 100.0e3);
    r rcomp("r2", n7, n23, 1.5e3);
    c cicomp("cicomp", n23, n2, 5.0e-9);
    c cocomp("cocomp", n18, n15, 0.2e-9);
    qnl q1("q1", n2, n30, n3);
    qnl q2("q2", n4, n31, n4);
    qnl q3("q3", n19, n6, n5);
    qnl q4("q4", n6, n4, n11);
    qnl q5("q5", n6, n11, n12);
    qnl q6("q6", n7, n13, n12);
    qnl q7("q7", n7, n2, n13);
    qnl q8("q8", n19, n7, n21);
    qnl q9("q9", n19, n17, n18);
    qnl q10("q10", n17, n15, n16);
    qnl q11("q11", n3, n8, n22);
    qnl q12("q12", n8, n8, n20);
    qnl q13("q13", n14, n14, n12);
    qpl q14("q14", n15, n12, n10);
    qpl q15("q15", n20, n17, n18);
    r r1("r1", n5, n2, 25.0e3);
    r r2("r2", n5, n4, 25.0e3);
    r r3("r3", n22, n20, 2.4e3);
    r r4("r4", n8, n9, 18.0e3);
    r r5("r5", n9, n12, 3.6e3);
    r r6("r6", n11, n14, 3.0e3);
    r r7("r7", n19, n6, 10.0e3);
    r r8("r8", n19, n7, 10.0e3);
    r r9("r9", n9, n10, 10.0e3);
    r r10("r10", n10, n18, 30.0e3);
    r r11("r11", n19, n17, 20.0e3);
    r r12("r12", n15, n16, 10.0e3);
    r r13("r13", n16, n20, 75.0);
    r r14("r14", n13, n14, 3.0e3);
    r r15("r15", n21, n10, 3.0e3);

    op::monitor = 1;
    op();

    //plot::output = *argv;
    //plot::add("n16");
    //dc::monitor = 1;
    //dc(vin, -0.2, 0.2, 0.001);
}
