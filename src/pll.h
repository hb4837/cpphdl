
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

#include "hdl.h"

MODEL(pd) {
private:
    node& in1;
    node& in2;
    node& out;
    node& gnd;
    parameter<double> k;
    quantity<across> vin1;
    quantity<across> vin2;
    quantity<across> vout;
    quantity<through> iin1;
    quantity<through> iin2;
    quantity<through> iout;
public:
    pd(string s, node& IN1, node& IN2, node& OUT, node& GND, const double K) 
    : model(s), 
      in1(IN1), 
      in2(IN2), 
      out(OUT), 
      gnd(GND), 
      k("gain", K), 
      vin1("vin1", in1, gnd), 
      vin2("vin2", in2, gnd), 
      vout("vout", out, gnd), 
      iin1("iin1", in1, gnd), 
      iin2("iin2", in2, gnd), 
      iout("iout", out, gnd)
    { }

    void parameters() { }

    void behaviour() {
        iin1 == 0;
        iin2 == 0;
        vout == k * vin1 * vin2;
    }
};

MODEL(lf) {
private:
    node& in;
    node& out;
    node& gnd;
    parameter<double> fp;
    quantity<across> vin;
    quantity<across> vout;
    quantity<through> iin;
    quantity<through> iout;
public:
    lf(string s, node& IN, node& OUT, node& GND, const double FP)
    : model(s), in(IN), out(OUT), gnd(GND), fp("fp", FP),
      vin("vin", in, gnd),
      vout("vout", out, gnd),
      iin("vin", in, gnd),
      iout("iout", out, gnd)
    { }

    void parameters() { }

    void behaviour() { 
        iin == 0.0;
        vin == vout + vout.dot() / fp;
    }
};

MODEL(vco) {
private:
    node& inp;
    node& inm;
    node& p;
    node& m;
    parameter<double> ampl;
    parameter<double> f0;
    parameter<double> gain;
    quantity<across> vin;
    quantity<across> vpm;
    quantity<through> ipm;
    quantity<double> integ;
    double dt;
public:
    vco(
        string s, 
        node& INP, node& INM, node& P, node& M, 
        const double AMPL, const double F0, const double GAIN
    ) : model(s), inp(INP), inm(INM), p(P), m(M), 
        ampl("amplitude", AMPL), 
        f0("center frequency", F0), 
        gain("frequency gain", GAIN), 
        vin("vin", inp, inm), 
        vpm("vpm", p, m), 
        ipm("ipm", p, m), 
        integ("integ")
    { 
        dt = 0.0;
    }

    void parameters() {
    }

    void behaviour() {
        quantity<double> ts;

        if(dc_domain()) {
            integ == 0.0;
        }
        else {
            vin == integ.dot();
            ts = 1.0 / (f0 + gain * vin) / 20.0;
            if(ts > 0.0)
                timestep(ts);
        }
        vpm == ampl * sin(constant::twopi * (f0 * time() + gain * integ));
    }
};

