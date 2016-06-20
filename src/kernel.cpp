
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

#include "global.h"
#include "kernel.h"
#include "node.h"
#include "quantity.h"
#include "equation.h"
#include "model.h"
#include "parameter.h"
#include "matrix.h"

tabulator indent;

string plot::output = "noname";

void plot::add(const string& s) {
    kernel::plot(s);
}

void plot::add(const string& s1, const string& s2) {
    kernel::plot(s1);
    kernel::plot(s2);
}

void plot::add(const string& s1, const string& s2, const string& s3) {
    kernel::plot(s1);
    kernel::plot(s2);
    kernel::plot(s3);
}

void plot::add(const string& s1, const string& s2, const string& s3, 
    const string& s4) 
{
    kernel::plot(s1);
    kernel::plot(s2);
    kernel::plot(s3);
    kernel::plot(s4);
}

void plot::add(const string& s1, const string& s2, const string& s3, 
    const string& s4, const string& s5) 
{
    kernel::plot(s1);
    kernel::plot(s2);
    kernel::plot(s3);
    kernel::plot(s4);
    kernel::plot(s5);
}

int op::monitor = 0;
int op::newton::itlim = 100;
double op::newton::abstol = 1.0e-12;
double op::newton::reltol = 1.0e-3;
double op::newton::vntol = 1.0e-6;

int op::ramp::src::itlim = 0;
double op::ramp::src::step = 1.0e-3;

int dc::monitor = 0;

int tran::itlim = 10;
int tran::monitor = 0;
int tran::method = trapezoidal;
int tran::order = 2;
int tran::stepsize = step_variable;
double tran::trtol = 7.0;
double tran::chgtol = 1.0e-14;
double tran::tsmin = 1.0e-15;
double tran::tscale = 1.0e6;

list<node *> kernel::nodes = (list<node *>) 0;
list<model *> kernel::models = (list<model *>) 0;
//list<equation *> kernel::eqs = (list<equation *>) 0;
list<quantity_ *> kernel::quantities = (list<quantity_ *>) 0;
model *kernel::curModel = 0;
int kernel::iterno = 0;
int kernel::analysis = NONE;
int kernel::status = NONE;
bool kernel::conv = false;
double kernel::stepFactor = 1.0;
matrix kernel::g = matrix();
//matrix kernel::a = matrix();
string kernel::plFile = "";
double kernel::delta[3] = {0.0, 0.0, 0.0};
double kernel::t = 0.0;
int kernel::order = 0;
list<double> kernel::breakpoints = (list<double>) 0;
bool kernel::setupFlag = false;
bool kernel::limited = false;
double kernel::userTimeStep = 1.0e9;
bool kernel::icFlag = false;
double kernel::f = 0.0;
int kernel::numNodes = 0;
int kernel::numQuantities = 0;
int kernel::nerrors = 0;
int kernel::nwarnings = 0;
int kernel::nfatal = 0;

void kernel::insert(node *n) {
    nodes.push_back(n);
    //n->index = nodes.size() - 1;
}

void kernel::insert(model *m) {
    models.push_back(m);
    curModel = m;
}

void kernel::insert(quantity_ *q) {
    quantities.push_back(q);
    //q->index = nodes.size() + quantities.size() - 1;
}

//void kernel::insert(equation *e) {
    //eqs.push_back(e);
    //e->index = nodes.size() + eqs.size() - 1;
//}

model *kernel::getCurModel() {
    return(curModel);
}

void kernel::evaluate() {
    list<model *>::iterator im;

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        (*im)->evaluate();
    }
}

void kernel::acLoad() {
}

void kernel::init() {
    // list<node *>::iterator in;
    // list<model *>::iterator im;

    status = INIT;

    // Calculate the node indices
    for(auto in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        (*in)->init();
    }

    for(auto im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        (*im)->init();
    }

    status = NONE;

    if(nerrors > 0) {
        exit(1);
    }
}

void kernel::load() {
    list<model *>::iterator im;

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        (*im)->load();
    }
}

void kernel::parameters() {
    list<model *>::iterator im;
    list<equation>::iterator ieq;
    //list<parameter_ *>::iterator ip;
    //list<quantity_ *>::iterator iq;

    //cerr << "----------------" << endl;
    //cerr << "Model parameters" << endl;
    //cerr << "----------------" << endl;
    //cerr << endl;

    for(im = models.begin(); im != models.end(); ++im) {
        (*im)->parameters();
        if(!setupFlag) {
            for(ieq = (*im)->equations.begin(); ieq != (*im)->equations.end();
                ++ieq)
            {
                (*ieq).setup();
            }
        }

        //cerr << (*im)->name << ":" << endl;
        //for(unsigned int i = 0; i <= (*im)->name.length(); ++i)
            //cerr << "-";
        //cerr << endl;
        //for(ip = (*im)->params.begin(); ip != (*im)->params.end(); ++ip) {
            //cerr.setf(ios::left);
            //cerr << setfill('.');
            //cerr << setw(10) << ((*ip)->name).data();
            //cerr.unsetf(ios::left);
            //cerr << setfill(' ');
            //cerr << setw(15) << (*ip)->getValueString();
            //cerr << endl;
        //}
        //cerr << endl;
    }

    setupFlag = true;
}

void kernel::op() {
    list<node *>::iterator in;
    list<quantity_ *>::iterator iq;

    parameters();
    init();

    iterno = 1;

    cerr << "---------------" << endl;
    cerr << "Operating point" << endl;
    cerr << "---------------" << endl;
    cerr << endl;

    iterno = iterate(op::newton::itlim, op::monitor);
    if(icFlag) {
        icFlag = false;
        iterno = iterate(op::newton::itlim, op::monitor);
    }

    if(!conv && op::ramp::src::itlim != 0) {
        cerr << "*op*: newton algorithm failed to converge" << endl;
        cerr << "      trying source stepping..." << endl;
        iterno = sourceStepping(op::ramp::src::itlim);
        if(conv && iterno <= op::ramp::src::itlim)
            cerr << "*message*: source stepping algorithm succeeded." << endl;
        else
            error("source stepping algorithm failed.");
    }

    if(!conv) {
        cerr << "*op*: no convergence in op analysis" << endl;
        cerr << "      last solution vector:" << endl;
        for(in = nodes.begin(); in != nodes.end(); ++in) {
            if((*in)->ground)
                continue;
            cerr << "         " << ((*in)->name).data() << " = ";
            cerr << setw(15) << (*in)->value << endl;
        }
        for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
            cerr << "         " << ((*iq)->name).data() << " = ";
            cerr << setw(15) << (*iq)->value << endl;
        }
        exit(1);
    }

    if(analysis == OP) {
        cerr << endl << "*op*: operating point analysis results:" << endl;
        for(in = nodes.begin(); in != nodes.end(); ++in) {
            if((*in)->ground)
                continue;
            cerr << setw(15) << ((*in)->name).data();
            cerr << setw(15) << (*in)->value << endl;
        }
        //for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
            //cerr << setw(15) << ((*iq)->name).data();
            //cerr << setw(15) << (*iq)->value << endl;
        //}
    }
    status = OP_DONE;

    cerr << endl;
    cerr << "Operating point analysis finished" << endl;
    cerr << endl;
    cerr << "Statistics: " << endl;
    cerr << "    number of iterations:   " << setw(15) << iterno << endl;
    cerr << "    number of unknowns:     " << setw(15) << g.getSize() << endl;
    cerr << endl;

}

void op() {
    kernel::analysis = kernel::OP;
    kernel::op();
    kernel::analysis = kernel::NONE;
}

void kernel::dc(model& m, const double& swb, const double& swe, 
    const double& sws)
{
    list<node *>::iterator in;
    list<quantity_ *>::iterator iq;
    int numdc = 0;
    int numnit = 0;
    double percent;

    plFile = plot::output + ".dc";
    ofstream pl(plFile.c_str());

    parameters();
    init();

    cerr << "--------------------" << endl;
    cerr << "DC transfer analysis" << endl;
    cerr << "--------------------" << endl;
    cerr << endl;
    if(dc::monitor > 0) {
        cerr << setw(15) << "sweep value";
        cerr << setw(5) << "it";
        cerr << setw(8) << "percent";
        cerr << endl;
    }

    pl << "# sweep";
    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        if((*in)->printFlag)
            pl << " " << ((*in)->name).data();
    }
    for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
        if((*iq)->printFlag)
            pl << " " << ((*iq)->name).data();
    }
    pl << endl;

    double fac = 100.0 / (swe - swb);

    for(double val = swb; val <= swe; val += sws) {
        percent = fac * (val - swb);
        numdc++;

        m.sweep(val);
        iterno = 1;

        iterno = iterate(op::newton::itlim, 0);
        numnit += iterno;
        //if(icFlag) {
            //icFlag = false;
            //iterno = iterate(op::newton::itlim, dc::monitor);
        //}

        if(dc::monitor > 0) {
            if(numdc % dc::monitor == 0) {
                cerr << setw(15) << val;
                cerr << setw(5) << iterno;
                cerr << setprecision(3) << setw(8) << percent;
                cerr << endl;
            }
        }

        if(!conv) {
            cerr << "*op*: no convergence in dc analysis" << endl;
            cerr << "      last sweep value:" << val << endl;
            cerr << "      last solution vector:" << endl;
            for(in = nodes.begin(); in != nodes.end(); ++in) {
                if((*in)->ground)
                    continue;
                cerr << "         " << ((*in)->name).data() << " = ";
                cerr << setw(15) << (*in)->value << endl;
            }
            for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
                cerr << "         " << ((*iq)->name).data() << " = ";
                cerr << setw(15) << (*iq)->value << endl;
            }
            exit(1);
        }

        pl.setf(ios::scientific);
        pl << val;
        for(in = nodes.begin(); in != nodes.end(); ++in) {
            if((*in)->ground)
                continue;
            if((*in)->printFlag)
                pl << " " << (*in)->value;
        }
        for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
            if((*iq)->printFlag)
                pl << " " << (*iq)->value;
        }
        pl << endl;
    }

    cerr << endl;
    cerr << "DC transfer analysis finished" << endl;
    cerr << endl;
    cerr << "Statistics: " << endl;
    cerr << "    total number of sweep steps: " << setw(15) << numdc << endl;
    cerr << "    total number of iterations:  " << setw(15) << numnit << endl;
    cerr << endl;

}

void kernel::tran(const double& ts, const double& te, const double& tb,
    const double& tsmax)
{
    int i;
    iterno = 0;
    int numnit = 0;
    int numtp = 0;
    int numrej = 0;
    double delmin;
    double delnew;
    list<node *>::iterator in;
    list<quantity_ *>::iterator iq;
    double minTs = te;
    double tfac = 100.0 / (te - tb);
    double percent;
    double lte = 0.0;

    plFile = plot::output + ".tran";
    ofstream pl(plFile.c_str());

    pl << "# time";
    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        if((*in)->printFlag)
            pl << " " << ((*in)->name).data();
    }
    for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
        if((*iq)->printFlag)
            pl << " " << ((*iq)->name).data();
    }
    pl << endl;

    delta[0] = dmin(ts, 0.02 * (te - tb));
    delmin = tran::tsmin;

    if(status != OP_DONE) {
        op();
    }

    acceptTimeStep();

#ifdef DEBUG
            cout << endl << endl << endl;
            cout << "===============================================" << endl;
            cout << "Starting transient analysis" << endl;
            cout << "===============================================" << endl;
#endif

    kernel::analysis = TRAN;

    cerr << endl;
    cerr << "------------------" << endl;
    cerr << "Transient analysis" << endl;
    cerr << "------------------" << endl;
    cerr << endl;
    if(tran::monitor > 0) {
        cerr << setw(5) << "+/-";
        cerr << setw(15) << "time";
        cerr << setw(5) << "it";
        cerr << setw(15) << "tstep";
        cerr << setw(15) << "lte";
        cerr << setw(6) << "order";
        cerr << setw(8) << "percent";
        cerr << endl;
    }

    for(i = 1; i < 3; ++i)
        delta[i] = delta[0];

    t = delta[0];
    percent = t * tfac;
    order = 1;

    do {
        do {
#ifdef DEBUG
            cout << endl << endl << endl;
            cout << "=====================================" << endl;
            cout << "Time: " << setw(15) << t << endl;
            cout << "=====================================" << endl;
#endif
            if(delta[0] < delmin) {
                error("internal timestep too small in transient analysis");
                exit(1);
            }

            iterno = iterate(tran::itlim, 0);
            numnit += iterno;

            // no convergence in iterate
            if(!conv) {
                if(tran::monitor > 0) {
                    pl.setf(ios::scientific);
                    cerr << setw(5) << "?";
                    cerr << setprecision(4) << setw(15) << t;
                    cerr << setw(5) << iterno;
                    cerr << setprecision(4) << setw(15) << delta[0];
                    cerr << setw(15) << "";
                    cerr << setw(6) << order;
                    cerr << setprecision(3) << setw(8) << percent;
                    cerr << endl;
                }

                t -= delta[0];
                delta[0] *= 0.125;
                t += delta[0];
                percent = t * tfac;
                order = 1;
                numrej++;

#ifdef DEBUG
                cout << "iterate() did not converge" << endl;
                cout << "trying again with tstep = " << delta[0];
#endif

                rejectTimeStep();
            }
        } while(!conv);

        delnew = delta[0];
        if(numtp > 0)
            lte = trunc(delnew);
        else
            delnew = delta[0];

        if(delnew >= 0.9 * delta[0]) {
            if(tran::monitor > 0) {
                if(numtp % tran::monitor == 0) {
                    pl.setf(ios::scientific);
                    cerr << setw(5) << "+";
                    cerr << setprecision(4) << setw(15) << t;
                    cerr << setw(5) << iterno;
                    cerr << setprecision(4) << setw(15) << delta[0];
                    cerr << setprecision(4) << setw(15) << lte;
                    cerr << setw(6) << order;
                    cerr << setprecision(3) << setw(8) << percent;
                    cerr << endl;
                }
            }
            order = tran::order;
            numtp++;

            delnew = dmin(delnew, 2.0 * delta[0]);
            delnew = dmin(delnew, tsmax);
            delnew = dmin(delnew, te - t);
            delnew = dmin(delnew, userTimeStep);

            delta[2] = delta[1];
            delta[1] = delta[0];
            delta[0] = delnew;

#ifdef DEBUG
            cout << "Setting tstep = " << delnew;
            cout << endl << "Accepting time step" << endl;
#endif
            acceptTimeStep();

            // write plot data
            pl.setf(ios::scientific);
            pl << t * tran::tscale;
            for(in = nodes.begin(); in != nodes.end(); ++in) {
                if((*in)->ground)
                    continue;
                if((*in)->printFlag)
                    pl << " " << (*in)->value;
            }
            for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
                if((*iq)->printFlag)
                    pl << " " << (*iq)->value;
            }
            pl << endl;

            t += delta[0];
            percent = t * tfac;
        }
        else {
            // truncation error too large
            if(tran::monitor > 0) {
                if(numtp % tran::monitor == 0) {
                    pl.setf(ios::scientific);
                    cerr << setw(5) << "-";
                    cerr << setprecision(4) << setw(15) << t;
                    cerr << setw(5) << iterno;
                    cerr << setprecision(4) << setw(15) << delta[0];
                    cerr << setprecision(4) << setw(15) << lte;
                    cerr << setw(6) << order;
                    cerr << setprecision(3) << setw(8) << percent;
                    cerr << endl;
                }
            }
            order = tran::order;

            t -= delta[0] - delnew;
            percent = t * tfac;
            delta[0] = delnew;
            numrej++;
        }

        minTs = dmin(minTs, delta[0]);
    } while(t < te);

    cerr << endl;
    cerr << "Transient analysis finished" << endl;
    cerr << endl;
    cerr << "Statistics: " << endl;
    cerr << "    total number of timesteps:   " << setw(15) << numtp << endl;
    cerr << "    number of rejected timesteps:" << setw(15) << numrej << endl;
    cerr << "    total number of iterations:  " << setw(15) << numnit << endl;
    cerr << "    minimum timestep:            " << setw(15) << minTs << endl;
    cerr << endl;

    status = TRAN_DONE;
    kernel::analysis = NONE;
}

double kernel::trunc(double& delnew) {
    double x = 0.0;
    double a0 = 0.0, a2 = 0.0;
    double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
    double q0 = 0.0, q1 = 0.0, q2 = 0.0, q3 = 0.0;
    double i0 = 0.0, i1 = 0.0;
    double tol1 = 0.0, tol2 = 0.0;
    double diff = 0.0;
    double newStep = 0.0;
    list<model *>::iterator im;
    //list<equation>::iterator ieq;
    //map<const string, quantity_ *>::const_iterator iq;
    list<quantity_ *>::iterator iq;

    if(tran::stepsize == step_fixed) {
        delnew = 1.0e-9;
        return(0.0);
    }

    double delold = 2.0 * delnew;

    if(tran::method == trapezoidal) {
        switch(order) {
            case 1:
                x = delta[0] + delta[1];
                d0 = 1.0 / (delta[0] * x);
                d2 = 1.0 / (delta[1] * x);
                d1 = -d0 - d2;
                break;
            case 2:
                x = delta[0] + delta[1] + delta[2];
                a0 = 1.0 / ((delta[0] + delta[1]) * x);
                a2 = 1.0 / ((delta[1] + delta[2]) * x);
                d0 = a0 / delta[0];
                d3 = -a2 / delta[2];
                x = (-a0 - a2) / delta[1];
                d1 = x - d0;
                d2 = -d3 - x;
                break;
        }
    }

    /*for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        for(ieq = (*im)->equations.begin(); ieq != (*im)->equations.end();
            ++ieq) 
        {
            for(iq = (*ieq).quantities.begin(); iq != (*ieq).quantities.end();
                ++iq) 
            {
                if(!iq->second->dynamicFlag)
                    continue;

                q0 = iq->second->trValue[0];
                q1 = iq->second->trValue[1];
                q2 = iq->second->trValue[2];
                if(order == 2)
                    q3 = iq->second->trValue[3];

                diff = d0 * q0 + d1 * q1 + d2 * q2 + op::newton::abstol;
                if(order == 2)
                    diff += d3 * q3;

                if(diff != 0.0) {
                    i0 = iq->second->dotValue;
                    i1 = iq->second->trDotValue[0];
                    tol1 = op::newton::reltol * dmax(fabs(i0), fabs(i1)) + 
                        op::newton::abstol;
                    tol2 = op::newton::reltol * (dmax(fabs(q0), fabs(q1)) +
                        tran::chgtol) / delta[0];
                    newStep = tran::trtol * dmax(tol1, tol2) / fabs(diff);

                    if(order == 2)
                        newStep = sqrt(2.0 * newStep);
                    delnew = dmin(newStep, delold);
                }
            }
        }
    }*/
    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        for(iq = (*im)->quantities.begin(); iq != (*im)->quantities.end();
            ++iq) 
        {
            //if(!(*iq)->dynamicFlag)
                //continue;

            q0 = (*iq)->trValue[0];
            q1 = (*iq)->trValue[1];
            q2 = (*iq)->trValue[2];
            if(order == 2)
                q3 = (*iq)->trValue[3];

            diff = d0 * q0 + d1 * q1 + d2 * q2 + op::newton::abstol;
            if(order == 2)
                diff += d3 * q3;

            if(diff != 0.0) {
                i0 = (*iq)->dotValue;
                i1 = (*iq)->trDotValue[0];
                tol1 = op::newton::reltol * dmax(fabs(i0), fabs(i1)) + 
                    op::newton::abstol;
                tol2 = op::newton::reltol * (dmax(fabs(q0), fabs(q1)) +
                    tran::chgtol) / delta[0];
                newStep = tran::trtol * dmax(tol1, tol2) / fabs(diff);

                //cerr << "!" << (*iq)->name << " " << newStep << endl;

                if(order == 2)
                    newStep = sqrt(2.0 * newStep);
                delnew = dmin(newStep, delold);
            }
        }
    }
    return(0.0);
}

void kernel::acceptTimeStep() {
    list<node *>::iterator in;
    list<model *>::iterator im;
    //list<equation>::iterator ieq;
    //map<const string, quantity_ *>::const_iterator miq;
    list<quantity_ *>::iterator liq;

    for(in = nodes.begin(); in != nodes.end(); ++in) {
        (*in)->trValue[2] = (*in)->trValue[1];
        (*in)->trValue[1] = (*in)->trValue[0];
        (*in)->trValue[0] = (*in)->value;
    }

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
//        for(ieq = (*im)->equations.begin(); ieq != (*im)->equations.end();
//            ++ieq) 
//        {
//            for(miq = (*ieq).quantities.begin(); miq != (*ieq).quantities.end();
//                ++miq) 
//            {
//                if(order == 2)
//                    miq->second->trValue[3] = miq->second->trValue[2];
//                miq->second->trValue[2] = miq->second->trValue[1];
//                miq->second->trValue[1] = miq->second->trValue[0];
//                miq->second->trValue[0] = miq->second->value;
//                if(miq->second->dynamicFlag) {
//                    miq->second->trDotValue[1] = miq->second->trDotValue[0];
//                    miq->second->trDotValue[0] = miq->second->dotValue;
//                }
//            }
//        }
        for(liq = (*im)->quantities.begin(); liq != (*im)->quantities.end(); 
            ++liq) 
        {
            if(order == 2)
                (*liq)->trValue[3] = (*liq)->trValue[2];
            (*liq)->trValue[2] = (*liq)->trValue[1];
            (*liq)->trValue[1] = (*liq)->trValue[0];
            (*liq)->trValue[0] = (*liq)->value;
            if((*liq)->dynamicFlag) {
                (*liq)->trDotValue[1] = (*liq)->trDotValue[0];
                (*liq)->trDotValue[0] = (*liq)->dotValue;
            }
            //if((*liq)->name == "q(c1)")
                //cout << "# " << **liq << endl;
        }
    }
}

void kernel::rejectTimeStep() {
    list<node *>::iterator in;
    list<model *>::iterator im;
    //list<equation>::iterator ieq;
    //map<const string, quantity_ *>::const_iterator miq;
    list<quantity_ *>::iterator liq;

    for(in = nodes.begin(); in != nodes.end(); ++in) {
        (*in)->value = (*in)->trValue[0];
    }

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
//        for(ieq = (*im)->equations.begin(); ieq != (*im)->equations.end();
//            ++ieq) 
//        {
//            for(miq = (*ieq).quantities.begin(); miq != (*ieq).quantities.end();
//                ++miq) 
//            {
//                miq->second->value = miq->second->trValue[0];
//                //miq->second->dotValue = miq->second->trDotValue[0];
//            }
//        }
        for(liq = (*im)->quantities.begin(); liq != (*im)->quantities.end();
            ++liq) 
        {
            (*liq)->value = (*liq)->trValue[0];
        }
    }
}

void discontinuity(const double& bp) {
    kernel::breakpoints.push_back(bp);
    kernel::breakpoints.sort();
}

bool dc_domain() {
    if(kernel::analysis == kernel::OP || kernel::analysis == kernel::DC) 
        return true;
    else
        return false;
}

bool ac_domain() {
    if(kernel::analysis == kernel::AC) 
        return true;
    else
        return false;
}

bool time_domain() {
    if(kernel::analysis == kernel::TRAN) 
        return true;
    else
        return false;
}

void timestep(const double& ts) {
    kernel::userTimeStep = dmin(kernel::userTimeStep, ts);
}

void timestep(const quantity_& q) {
    kernel::userTimeStep = dmin(kernel::userTimeStep, q.current());
}

void dc(model& m, const double& swb, const double& swe, const double& sws) {
    kernel::analysis = kernel::DC;
    kernel::dc(m, swb, swe, sws);
    kernel::analysis = kernel::NONE;
}

void kernel::ac(const acscale& acs, const int& n, const double& fb, 
    const double& fe) 
{
    double fs;

    g.clear();

    for(f = fb; f <= fe; f += fs) {
        acLoad();
        g.factor();
        g.solve();
        getSolution();
    }
}

void ac(const acscale& acs, const int& n, const double& fb, const double& fe) {
    kernel::analysis = kernel::AC;
    ac(acs, n, fb, fe);
    kernel::analysis = kernel::NONE;
}

void tran(const double& ts, const double& te) {
    double tsmax = dmin(0.02 * te, ts);
    kernel::tran(ts, te, 0.0, tsmax);
}

void tran(const double& ts, const double& te, const double& tb) {
    double tsmax = dmin(0.02 * (te - tb), ts);
    kernel::tran(ts, te, tb, tsmax);
}

void tran(const double& ts, const double& te, const double& tb, 
    const double& tsmax)
{
    kernel::tran(ts, te, tb, tsmax);
}

int kernel::iterate(const int& itlim, const int& mon) {
    int iterno = 0;

    if(mon > 0) {
        cerr << endl;
        cerr << setw(5) << "it";
        cerr << setw(5) << "+/-";
        cerr << endl;
    }

    do {
#ifdef DEBUG
        cout << endl << endl << endl;
        cout << "===========================" << endl;
        cout << "Starting iteration no." << setw(5) << iterno + 1 << endl;
        cout << "===========================" << endl;
#endif
        g.clear();

        evaluate();

        load();

#ifdef DEBUG
        g.info();
#endif
        g.factor();
        g.solve();

        getSolution();

        if(nfatal > 0)
            exit(1);

        conv = convergence();

        iterno++;

        if(mon > 0) {
            if(iterno % mon == 0) {
                cerr << setw(5) << iterno;
                if(conv)
                    cerr << setw(5) << "+";
                else {
                    if(limited) {
                        cerr << setw(5) << "*";
                        limited = false;
                    }
                    else {
                        cerr << setw(5) << "-";
                    }
                }
                cerr << endl;
            }
        }
        else
            limited = false;
    } while(!conv && iterno <= itlim);

    return(iterno);
}

int kernel::sourceStepping(const int& itlim) {
    int totalIterations = 0;
    double step = 1.0e-3;
    list<node *>::iterator in;
    list<quantity_ *>::iterator iq;
    list<model *>::iterator im;

    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        (*in)->reset();
    }

    for(im = models.begin(); im != models.end(); ++im) {
        for(iq = (*im)->quantities.begin(); iq != (*im)->quantities.end();
            ++iq) 
        {
            (*iq)->reset();
        }
    }

    iterno = 0;

    stepFactor = 0.0;
    iterno = iterate(op::newton::itlim, op::monitor);

    cerr << endl << "  factor iterations total conv" << endl;
    do {
        stepFactor += step;

        iterno = 0;
        iterno = iterate(op::newton::itlim, op::monitor);

        totalIterations += iterno;

        cerr << setw(8) << stepFactor;
        cerr << setw(11) << iterno;
        cerr << setw(6) << totalIterations;
        if(conv)
            cerr << setw(5) << "+" << endl;
        else {
            cerr << setw(5) << "-" << endl;
            stepFactor -= step;
        }

        if(!conv || iterno >= 6)
            step *= 0.5;
        else if(conv && (iterno <= 3))
            step *= 2.0;

        if(step + stepFactor >= 1.0)
            step = 1.0 - stepFactor;

    } while(stepFactor < 1.0 && totalIterations <= itlim);

    cerr << endl;
    return(totalIterations);
}

bool kernel::convergence() {
    list<node *>::iterator in;
    list<model *>::iterator im;

#ifdef DEBUG
    cout << endl;
    cout << "---------------------" << endl;
    cout << "kernel::convergence()" << endl;
    cout << "---------------------" << endl;
    cout << endl;
    //cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
    cout << "  node/quantity          error          limit" << endl;
    //cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
#endif

    if(limited) {
        return(false);
    }

    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        if(!(*in)->convergence())
            return(false);
    }

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        if(!(*im)->convergence())
            return(false);
    }

//#ifdef DEBUG
//    cout << setfill('-') << setw(60) << "" << setfill(' ') << endl;
//#endif
    return(true);
}

void kernel::getSolution() {
    list<node *>::iterator in;
    list<model *>::iterator im;
    list<quantity_ *>::iterator iq;

#ifdef DEBUG
    cout << endl;
    cout << "---------------------" << endl;
    cout << "kernel::getSolution()" << endl;
    cout << "---------------------" << endl << endl;
    cout << endl;
    cout << "  node/quantity          value          delta" << endl;
#endif

    // node values
    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        (*in)->getSolution();
    }

    // quantity values
    // first we need to calculate those quantities 
    // which are NOT part of the solution vector  (== explicit quantities)
    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        for(iq = (*im)->quantities.begin(); iq != (*im)->quantities.end(); 
            ++iq) 
        {
            (*iq)->getSolution();
        }
    }

    for(im = models.begin(); im != models.end(); ++im) {
        curModel = *im;
        (*im)->getSolution();
    }
}

void kernel::plot(const string& s) {
    list<node *>::iterator in;
    //list<model *>::iterator im;
    //list<equation>::iterator ieq;
    list<quantity_ *>::iterator liq;
    //map<const string, quantity_ *>::const_iterator iq;
    bool found = false;

    for(in = nodes.begin(); in != nodes.end(); ++in) {
        if((*in)->ground)
            continue;
        if((*in)->name == s) {
            (*in)->printFlag = true;
            found = true;
        }
    }

//    for(im = models.begin(); im != models.end(); ++im) {
//        for(ieq = (*im)->equations.begin(); ieq != (*im)->equations.end(); 
//            ++ieq)
//        {
//            if(!setupFlag)
//                (*ieq).setup();
//
//            for(iq = (*ieq).quantities.begin(); iq != (*ieq).quantities.end();
//                ++iq) 
//            {
//                if(iq->second->name == s) {
//                    iq->second->printFlag = true;
//                    found = true;
//                }
//            }
//        }
//    }
            for(liq = quantities.begin(); liq != quantities.end(); ++liq) {
                if((*liq)->name == s) {
                    (*liq)->printFlag = true;
                    found = true;
                }
            }
    setupFlag = true;

    if(!found)
        warning({"output variable ", s.c_str(), " does not exist"});
}

void plot(char *args, ...) {
    va_list ap;
    char *p = args;

    kernel::plot(p);

    va_start(ap, args);
    for(;;) {
        p = va_arg(ap, char *);
        if(p == 0) 
            break;
        kernel::plot(p);
    }
    va_end(ap);
}

void print(const int& n, ...) {
    va_list ap;
    int i;
    char *p;

    va_start(ap, n);
    for(i = 0; i < n; ++i) {
        p = va_arg(ap, char *);
        kernel::plot(p);
    }
    va_end(ap);
}

/*
void error(const string& s) {
    cerr << "*error*: ";
    cerr << s << endl;
    ++kernel::nerrors;
}

void error(const string& s1, const string& s2) {
    cerr << "*error*: ";
    cerr << s1;
    cerr << s2 << endl;
    ++kernel::nerrors;
}

void error(const string& s1, const string& s2, const string& s3) {
    cerr << "*error*: ";
    cerr << s1;
    cerr << s2;
    cerr << s3 << endl;
    ++kernel::nerrors;
}
*/

/*
void error(const char *fmt, ...) {
    // const char *fmt = format.c_str();    
    printf("ERROR: ");
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    ++kernel::nerrors;
}
*/

void error(const string& msg) {
    cerr << "ERROR: " << msg << endl;
}

void error(initializer_list<string> err) {
    cerr << "ERROR:";
    for (auto& s: err)
        cerr << " " << s;
    cerr << endl;
}

void warning(const string& msg) {
    cerr << "WARNING: " << msg << endl;
}

void warning(initializer_list<string> err) {
    // const char *fmt = format.c_str();    
    /*
    printf("WARNING: ");
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    */
    ++kernel::nwarnings;
}

void fatal(const string& msg) {
    cerr << "FATAL ERROR: " << msg << endl;
}

void fatal(initializer_list<string> err) {
    // const char *fmt = format.c_str();    
    /*
    printf("FATAL ERROR: ");
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    */
    ++kernel::nfatal;
}

