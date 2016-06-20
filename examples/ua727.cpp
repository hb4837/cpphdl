
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
    node n24("n24");
    node n25("n25");
    node n26("n26");
    node n27("n27");
    node n28("n28");
    node n29("n29");
    node n31("n31");
    node n32("n32");
    node n33("n33");
    node n34("n34");
    node n35("n35");
    node n36("n36");
    node n40("n40");

    vdc vcc1("vcc1", n34, gnd, 15.0);
    vdc vcc2("vcc1", n35, gnd, 15.0);
    vdc vee("vee", n36, gnd, -15.0);
    vdc vin("vin", n40, gnd, 0.0);

    r rs1("rs1", n40, n1, 1.0e3);
    r rs2("rs2", n12, gnd, 1.0e3);
    idc iz1("iz1", n36, n9, 0.62);
    r rz1("rz1", n36, n9, 10.0);
    idc iz2("iz2", n32, n31, 0.62);
    r rz2("rz2", n32, n31, 10.0);
    r r1("r1", n9, n31, 1.0e3);
    r r2("r2", n28, n9, 21.0e3);
    r r3("r3", n28, n19, 4.8e3);
    r r4("r4", n32, n36, 2.4e3);
    r r5("r5", n33, n36, 10.0);
    r r6("r6", n26, n19, 2.0e3);
    r r7("r7", n25, n36, 1.5e3);
    r r8("r8", n20, n36, 120.0e3);
    r r9("r9", n11, n3, 60.0e3);
    r r10("r10", n6, n8, 60.0e3);
    r r11("r11", n34, n5, 3.0e3);
    r r12("r12", n8, n9, 10.0e3);
    r r13("r13", n22, n36, 15.0e3);
    r r14("r14", n21, n36, 15.0e3);
    r r15("r15", n23, n36, 15.0e3);
    r r16("r16", n17, n9, 10.0e3);
    r r17("r17", n34, n15, 3.0e3);
    r r18("r18", n16, n17, 60.0e3);
    r r19("r19", n11, n14, 60.0e3);
    r r21("r21", n24, n36, 120.0e3);
    r rtemp("rtemp", n35, n29, 330.0e3);

    qnl q3("q3", n35, n29, n31);
    qnl q4("q4", n35, n32, n33);
    qnl q5("q5", n29, n33, n36);
    qnl q6("q6", n29, n28, n27);
    qnl q7("q7", n27, n27, n26);
    qnl q8("q8", n19, n19, n25);
    qnl q9("q9", n34, n1, n2);
    qnl q10("q10", n2, n19, n20);
    qnl q11("q11", n34, n34, n11);
    qnl q12("q12", n3, n2, n4);
    qnl q13("q13", n4, n19, n22);
    qpl q14("q14", n6, n3, n5);
    qnl q15("q15", n5, n6, n8);
    qnl q16("q16", n34, n8, n10);
    qnl q17("q17", n10, n19, n21);
    qnl q18("q18", n34, n17, n18);
    qnl q19("q19", n18, n19, n23);
    qnl q20("q20", n15, n16, n17);
    qpl q21("q21", n16, n14, n15);
    qnl q22("q22", n14, n13, n4);
    qnl q23("q23", n13, n19, n24);
    qnl q24("q24", n34, n12, n13);

    op::monitor = 1;
    op();
}
