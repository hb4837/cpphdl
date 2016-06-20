
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

int main(int argc, char *argv[]) {
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
    node n17("n17");
    node n18("n18");
    node n20("n20");
    node n21("n21");
    node n22("n22");
    node n23("n23");
    node n24("n24");
    node n25("n25");
    node n26("n26");
    node n27("n27");
    node n30("n30");

    vdc vcc("vcc", n27, gnd, 15.0);
    vdc vee("vee", n26, gnd, -15.0);
    //vdc vin("vin", n30, gnd, 0.0);
    vpulse vin("vin", n30, gnd, 0.0, 1.0, 1.0e-8, 1.0e-9);

    r rs1("rs1", n1, n30, 1.0e3);
    r rs2("rs2", n2, gnd, 1.0e3);
    r rs("rs", n24, n1, 100.0e3);
    r r1("r1", n10, n26, 1.0e3);
    r r2("r2", n9, n26, 50.0e3);
    r r3("r3", n11, n26, 1.0e3);
    r r4("r4", n12, n26, 3.0e3);
    r r5("r5", n15, n17, 39.0e3);
    r r6("r6", n21, n20, 40.0e3);
    r r7("r7", n14, n26, 50.0e3);
    r r8("r8", n18, n26, 50.0);
    r r9("r9", n24, n25, 25.0);
    r r10("r10", n23, n24, 50.0);
    r r11("r11", n13, n26, 50.0e3);
    c cc("cc", n22, n8, 30.0e-12);
    qnl q1("q1", n3, n2, n4);
    qnl q2("q2", n3, n1, n5);
    qpl q3("q3", n7, n6, n4);
    qpl q4("q4", n8, n6, n5);
    qnl q5("q5", n7, n9, n10);
    qnl q6("q6", n8, n9, n11);
    qnl q7("q7", n27, n7, n9);
    qnl q8("q8", n6, n15, n12);
    qnl q9("q9", n15, n15, n26);
    qpl q10("q10", n3, n3, n27);
    qpl q11("q11", n6, n3, n27);
    qpl q12("q12", n17, n17, n27);
    qnl q13("q13", n8, n13, n26);
    qpl q14("q14", n22, n17, n27);
    qnl q15("q15", n22, n22, n21);
    qnl q16("q16", n22, n21, n20);
    qnl q17("q17", n13, n13, n26);
    qnl q18("q18", n27, n8, n14);
    qnl q19("q19", n20, n14, n18);
    qnl q20("q20", n22, n23, n24);
    qpl q21("q21", n13, n25, n24);
    qnl q22("q22", n27, n22, n23);
    qpl q23("q23", n26, n20, n25);

    //op::monitor = 1;
    //op();

    plot::output = *argv;
    plot::add("n30", "n24");
    tran::order = 1;
    tran::monitor = 1;
    tran(1.0e-9, 1.0e-7);
}
