
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

class dmod : public d {
public:
    dmod(string s, node& P, node& M) : d(s, P, M) {
        rs = 40.0;
        tt = 0.1e-9;
        cjo = 0.9e-12;
    }
};

class qnd : public q {
public:
    qnd(string s, node& C, node& B, node& E) : q(s, C, B, E, "npn") {
        bf = 50.0;
        cje = 0.9e-12;
        cjc = 1.5e-12;
        vaf = 50.0;
        tf = 0.1e-9;
        tr = 10.0e-9;
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

    vdc vcc("vcc", n13, gnd, 5.0);
    vpulse vin("vin", n1, gnd, 0, 3.5, 1e-9, 1e-9, 1e-9, 40e-9);
    r rs("rs", n1, n2, 50);
    qnd q1("q1", n4, n3, n2);
    r rb1("rb1", n13, n3, 4e3);
    qnd q2("q2", n5, n4, n6);
    r rc2("rc2", n13, n5, 1.4e3);
    r re2("re2", n6, gnd, 1e3);
    qnd q3("q3", n7, n5, n8);
    r rc3("rc3", n13, n7, 100);
    dmod d1("d1", n8, n9);
    qnd q4("q4", n9, n6, gnd);
    qnd q5("q5", n11, n10, n9);
    r rb5("rb5", n13, n10, 4e3);
    dmod d2("d2", n11, n12);
    dmod d3("d3", n12, gnd);

    op::monitor = 1;
    op();

    plot::output = *argv;
    plot::add("n1", "n5", "n9");

    tran::order = 1;
    tran::monitor = 1;
    tran::itlim = 10;
    tran(1e-9, 100.0e-9);
}
