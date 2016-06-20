#include "hdl.h"

using namespace util;

class r : public model {
public:
    port pos;
    port neg;
    parameter<double> rnom;
    quantity<across> vpm;
    quantity<through> ipm;
    equation eq1;
public:

    r(const string s) : model(s), vpm("vpm", pos, neg), ipm("ipm", pos, neg), eq1(vpm, ipm) {
    }

    void parameters() {
        if (!rnom.defined())
            error({name, ": resistance value rnom must be defined"});
        if (rnom <= 0.0)
            error({name, ": non-positive resistance value"});
    }

    void load() {
        double vp = pos.getNode()->getValue();
        double vm = neg.getNode()->getValue();
        vpm.setValue(vp - vm);
        
        // ipm.setValue(vpm.getValue() / rnom);
        
        eq1.setF(ipm.getValue() - vpm.getValue() / rnom);
        eq1.setJ(vpm, -1.0 / rnom);
        eq1.setJ(ipm, 1.0);
    }
};

class c : public model {
public:
    port pos;
    port neg;
    parameter<double> value;
    quantity<across> vpm;
    quantity<through> ipm;
    quantity<double> q;

    c(string s) : model(s), vpm("vpm", pos, neg), ipm("ipm", pos, neg), q("q") {}

    void equations() {
        q == value * vpm;
        ipm == q.dot();
    }
};

class v : public model {
public:
    port pos;
    port neg;
    // DC
    parameter<double> dcval;
    // AC
    parameter<double> ampl;
    parameter<double> phase;
    // Pulse
    parameter<double> v1;
    parameter<double> v2;
    parameter<double> td;
    parameter<double> tr;
    parameter<double> tf;
    parameter<double> pw;
    parameter<double> per;
private:
    quantity<across> vpm;
    quantity<through> ipm;
    equation eq1;
private:
    int type;
    double mr, mf;
    double ttr;
    double ttw;
    double ttf;
    double ttper;
    double t;
    double val;
    static const int DC = 0;
    static const int AC = 1;
    static const int PULSE = 2;
    static const int PWL = 3;

public:

    v(string s) : model(s), vpm("vpm", pos, neg), ipm("ipm", pos, neg),
            eq1(vpm, ipm)
    {
    }

    void parameters() {
    }

    void dc(const double& v) {
        type = DC;
        dcval = v;
    }
    
    void pulse(const double& v1, const double& v2, const double& td,
            const double& tr, const double& tf, const double& pw=0.0, const double& per=0.0) {
        type = PULSE;
    }

    void load() {
        // f = vpm -val
        // df/dvpm = 1
        double vp = pos.getNode()->getValue();
        double vm = neg.getNode()->getValue();
        vpm.setValue(vp - vm);      
        
        eq1.setF(vpm.getValue() - dcval);
        eq1.setJ(vpm, 1.0);
    }
    
    void equations() {
        switch (type) {
            case DC:
                val = dcval;
                break;

            case AC:
                break;

            case PULSE:
                t = time();

                // reduce on first time period
                while (t > ttper)
                    t -= ttper;

                if (t > ttf)
                    val = v1; // after falling edge
                else if (t > ttw)
                    val = v1 + v2 - mf * (t - ttw); // during falling edge
                else if (t > ttr)
                    val = v1 + v2; // maximum value
                else if (t > td)
                    val = v1 + mr * (t - td); // during rising edge
                else
                    val = v1; // minimum value
                break;

            case PWL:
                break;

        }
        vpm == val;
    }

    void sweep(const double& d) {
        // ...
    }
};

