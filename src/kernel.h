
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

#ifndef _KERNEL_H
#define _KERNEL_H

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

using namespace std;

class tabulator {
private:
    string s;
public:
    string operator ++() { s += "  "; return(s); }
    string operator --() { s.erase(0,2); return(s); }

    friend ostream& operator <<(ostream& os, const tabulator& f) {
        os << f.s;
        return(os);
    }
};

/*
class constant {
public:
    static const double q = 1.6021918e-19;
    static const double epso = 8.854214871e-12;
    static const double epsox = 3.9;
    static const double epssi = 11.9;
    static const double k = 1.3806226e-23;
    static const double boltz = 1.3806226e-23;
    static const double ctok = 273.15;
    static const double koverq = 86.170869181e-6;

    static const double pi = 3.141592653589793240;
    static const double twopi = 6.283185307179586480;
    static const double sqrt2 = 1.41421356237;
    static const double e = 2.78182818285;
};
*/

class node;
class quantity_;
class model;
class equation;
class matrix;
class parameter_;

class plot {
public:
    static string output;
public:
    static void add(const string&);
    static void add(const string&, const string&);
    static void add(const string&, const string&, const string&);
    static void add(const string&, const string&, const string&, const string&);
    static void add(const string&, const string&, const string&, const string&,
        const string&);
};

class op {
public:
    static int monitor;

    class newton {
    public:
        static int itlim;
        static double abstol;
        static double reltol;
        static double vntol;
    };
    class ramp {
    public:
        class src {
        public:
            static int itlim;
            static double step;
        };
    };
};

class dc {
public:
    static int monitor;
};

static const int euler = 1;
static const int trapezoidal = 2;
static const int step_variable = 0;
static const int step_fixed = 1;

class tran {
public:
    static int itlim;
    static int method;
    static int order;
    static int monitor;
    static double trtol;
    static double chgtol;
    static double tsmin;
    static int stepsize;
    static double tscale;
};

enum acscale {DEC, OCT, LIN};

class kernel {
    friend class node;
    friend class model;
    friend class quantity_;
    friend class equation;
private:
    static list<node *> nodes;
    static list<quantity_ *> quantities;
    static list<model *> models;
    static model *curModel;
    static int analysis;
    static const int NONE = 0;
    static const int OP = 1;
    static const int DC = 2;
    static const int AC = 3;
    static const int TRAN = 4;
    static const int INIT = 5;
    static const int OP_DONE = 1;
    static const int TRAN_DONE = 2;
    static bool conv;
    static int iterno;
    static double stepFactor;
    static matrix g;
    static string plFile;
    static double t;
    static int order;
    static double delta[3];
    static list<double> breakpoints;
    static bool setupFlag;
    static double userTimeStep;
    static double f;
    static int numNodes;
    static int numQuantities;
    static int nerrors;
    static int nwarnings;
    static int nfatal;
public:
    static bool limited;
    static bool icFlag;
    static int status;
public:
    static void insert(node *);
    static void insert(quantity_ *);
    static void insert(model *);
    static model *getCurModel();
    static void evaluate();
    static void load();
    static void acLoad();
    static void init();
    static void parameters();
    static void op();
    static void dc(model&, const double&, const double&, const double&);
    static void ac(const acscale&, const int&, const double&, const double&);
    static void tran(const double&, const double&);
    static void tran(const double&, const double&, const double&);
    static void tran(const double&, const double&, const double&, const double&);
    static double trunc(double&);
    static void acceptTimeStep();
    static void rejectTimeStep();
    static int iterate(const int&, const int&);
    static int sourceStepping(const int&);
    static bool convergence();
    static void getSolution();
    static void plot(const string&);
    static double step() { return(stepFactor); }
    friend double step() { return(stepFactor); }
    friend double time() { return(t); }
    friend double freq() { return(f); }
    friend void op();
    friend void dc(model&, const double&, const double&, const double&);
    friend void ac(const acscale&, const int&, const double&, const double&);
    friend void tran(const double&, const double&, const double&, const double&);
    friend void plot(const string&, ...);
    friend void print(const int&, ...);
    friend void print(model&, quantity_&, ...);
    static double h(const int& i) { return delta[i]; }
    friend void discontinuity(const double& bt);
    static int currentOrder() { return order; }
    friend bool dc_domain();
    friend bool ac_domain();
    friend bool time_domain();
    friend void timestep(const double&);
    friend void timestep(const quantity_&);
    static bool initStatus() { return(status == INIT ? true : false); }
    friend void error(const string&);
    friend void error(initializer_list<string>);
    friend void warning(const string&);
    friend void warning(initializer_list<string>);
    friend void fatal(const string&);
    friend void fatal(initializer_list<string>);
};

#endif // _KERNEL_H
