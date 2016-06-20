#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <cmath>
#include <string>
// #include <iostream>
// #include <strstream>
#include "sim.h"

class quantity_;
class model;

using namespace std;

extern double gmin;
extern double temp;
extern double tempk;

#define dmin(a, b) (a < b ? a : b)
#define dmax(a, b) (a > b ? a : b)

// Friend functions from class kernel
//extern bool dc_domain();
//extern bool ac_domain();
//extern bool time_domain();
//extern void timestep(const double&);
//extern void timestep(const quantity_&);
//
//extern double step();
//extern double freq();
//extern void discontinuity(const double&);
//
//extern void op();
//extern void dc(model&, const double&, const double&, const double&);
//extern void tran(const double&, const double&);
//extern void tran(const double&, const double&, const double&);
//extern void tran(const double&, const double&, const double&, const double&);

extern double time();

#endif // GLOBAL_H
