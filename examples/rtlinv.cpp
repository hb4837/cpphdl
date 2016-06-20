
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

class qnd : public q {
public:
    qnd(string s, node& C, node& B, node& E) : q(s, C, B, E, "npn") {
        bf = 50.0;
        cje = 0.9e-12;
        cjc = 1.5e-12;
        vaf = 50.0;
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

    vdc vcc("vcc", n6, gnd, 5.0);
    vpulse vin("vin", n1, gnd, 0, 5, 2e-9, 2e-9, 2e-9, 10e-9);
    r rb1("rb1", n1, n2, 10e3);
    r rc1("rc1", n6, n3, 1e3);
    qnd q1("q1", n3, n2, gnd);
    r rb2("rb2", n3, n4, 10e3);
    qnd q2("q2", n5, n4, gnd);
    r rc2("rc2", n6, n5, 1e3);

    op::monitor = 1;
    op();

    plot::output = *argv;
    plot::add("n1", "n3", "n5");

    tran::tsmin = 1.0e-30;
    tran::monitor = 1;
    tran(0.2e-9, 20.0e-9);
}
