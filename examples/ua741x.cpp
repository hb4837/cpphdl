
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

class ua741 {
public:
    node n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, 
         n17, n18, n20, n21, n22, n23, n25;
    r rs1, rs2, rf, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
    c cc;
    qnl q1, q2;
    qpl q3, q4;
    qnl q5, q6, q7, q8, q9;
    qpl q10, q11, q12;
    qnl q13;
    qpl q14;
    qnl q15, q16, q17, q18, q19, q20;
    qpl q21;
    qnl q22;
    qpl q23;

    ua741(node& n30, node& n31, node& n24, node& n27, node& n26) :
        n1("n1"),
        n2("n2"),
        n3("n3"),
        n4("n4"),
        n5("n5"),
        n6("n6"),
        n7("n7"),
        n8("n8"),
        n9("n9"),
        n10("n10"),
        n11("n11"),
        n12("n12"),
        n13("n13"),
        n14("n14"),
        n15("n15"),
        n17("n17"),
        n18("n18"),
        n20("n20"),
        n21("n21"),
        n22("n22"),
        n23("n23"),
        n25("n25"),

        rs1("rs1", n1, n30, 1.0e3),
        rs2("rs2", n2, n31, 1.0e3),
        rf("rf", n24, n1, 100.0e3),
        r1("r1", n10, n26, 1.0e3),
        r2("r2", n9, n26, 50.0e3),
        r3("r3", n11, n26, 1.0e3),
        r4("r4", n12, n26, 3.0e3),
        r5("r5", n15, n17, 39.0e3),
        r6("r6", n21, n20, 40.0e3),
        r7("r7", n14, n26, 50.0e3),
        r8("r8", n18, n26, 50.0),
        r9("r9", n24, n25, 25.0),
        r10("r10", n23, n24, 50.0),
        r11("r11", n13, n26, 50.0e3),
        cc("cc", n22, n8, 30.0e-12),
        q1("q1", n3, n2, n4),
        q2("q2", n3, n1, n5),
        q3("q3", n7, n6, n4),
        q4("q4", n8, n6, n5),
        q5("q5", n7, n9, n10),
        q6("q6", n8, n9, n11),
        q7("q7", n27, n7, n9),
        q8("q8", n6, n15, n12),
        q9("q9", n15, n15, n26),
        q10("q10", n3, n3, n27),
        q11("q11", n6, n3, n27),
        q12("q12", n17, n17, n27),
        q13("q13", n8, n13, n26),
        q14("q14", n22, n17, n27),
        q15("q15", n22, n22, n21),
        q16("q16", n22, n21, n20),
        q17("q17", n13, n13, n26),
        q18("q18", n27, n8, n14),
        q19("q19", n20, n14, n18),
        q20("q20", n22, n23, n24),
        q21("q21", n13, n25, n24),
        q22("q22", n27, n22, n23),
        q23("q23", n26, n20, n25)
    { }
};

main(int argc, char *argv[]) {
    node gnd;
    node nvcc("nvcc");
    node nvee("nvee");
    node in("in");
    node out("out");

    vdc vcc("vcc", nvcc, gnd, 15.0);
    vdc vee("vee", nvee, gnd, -15.0);
    vdc vin("vin", in, gnd, 0.0);
    //vpulse vin("vin", in, gnd, 0.0, 1.0, 100.0e-6, 1.0e-9);
    ua741 ua741_1(in, gnd, out, nvcc, nvee);
    //c cl("cl", out, gnd, 1.0e-12);

    op::monitor = 1;
    op::ramp::src::itlim = 1000;
    op();

    //plot::output = *argv;
    //plot::add("in", "out");

    //dc::monitor = 1;
    //dc(vin, -0.2, 0.2, 0.002);
    //dc(vin, -15, 15, 0.1);

    //tran::monitor = 100;
    //tran(1.0e-6, 200.0e-6);
}
