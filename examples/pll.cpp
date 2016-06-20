
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
#include "pll.h"

class vpwl1 : public vpwl {
public:
    vpwl1(string s, node& P, node& M) : vpwl(s, P, M) {
        t[0]  = 0.0;        v[0] =  0.0;
//      t[0]  = 0.0;        v[0] = -4.0;
//      t[1]  = 10.0e-6;    v[1] = -4.0; 
//      t[2]  = 30.000e-6;  v[2] = -4.0; 
//      t[3]  = 30.001e-6;  v[3] = -3.0; 
//      t[4]  = 60.000e-6;  v[4] = -3.0; 
//      t[5]  = 60.001e-6;  v[5] = -2.0; 
//      t[6]  = 90.000e-6;  v[6] = -2.0; 
//      t[7]  = 90.001e-6;  v[7] = -1.0; 
//      t[8]  = 120.000e-6; v[8] = -1.0; 
//      t[9]  = 120.001e-6; v[9] =  0.0; 
//      t[10] = 150.000e-6; v[10] = 0.0; 
//      t[11] = 150.001e-6; v[11] = 1.0; 
//      t[12] = 180.000e-6; v[12] = 1.0; 
//      t[13] = 180.001e-6; v[13] = 2.0; 
//      t[14] = 210.000e-6; v[14] = 2.0; 
//      t[15] = 210.001e-6; v[15] = 3.0; 
//      t[16] = 240.000e-6; v[16] = 3.0; 
//      t[17] = 240.001e-6; v[17] = 4.0;
    }
};

int main(int argc, char *argv[]) {
    NODE(gnd);
    NODE(in);
    NODE(pll_in);
    NODE(pd_out);
    NODE(lf_out);
    NODE(vco_out);

    vpwl1 v1("v1", in, gnd);
    vco vco2("vco2", in, gnd, pll_in, gnd, 1.0, 7.0e6, 30.0e3);

    pd pd1("pd1", pll_in, vco_out, pd_out, gnd, 3.72);
    lf lf1("lf1", pd_out, lf_out, gnd, 700.0e3);
    vco vco1("vco1", lf_out, gnd, vco_out, gnd, 1.0, 7.0e6, 30.0e3);

    op::monitor = 1;
    op();

    plot::output = *argv;
    plot::add("in", "pll_in", "lf_out", "vco_out");
  
    tran::monitor = 1;
    tran::tsmin = 1.0e-21;
    tran::itlim = 100;
    tran(1.0e-7, 300.0e-6, 0.0, 1.0e-9);
}
