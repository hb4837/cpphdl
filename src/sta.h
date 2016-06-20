#ifndef STA_H
#define STA_H

#include <list>
#include <map>
#include <algorithm>
#include "sim.h"
#include "model.h"
#include "equation.h"
#include "node.h"
#include "quantity.h"
#include "matrix.h"

using namespace std;

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
//  class ramp {
//  public:
//      class src {
//      public:
//          static int itlim;
//          static double step;
//      };
//  };
};

class sta : public sim {
private:
    map<string, node *> nodes;
    map<string, quantity_ *> currents;
    double t;
    int numNodes;
    int numUnknowns;
    matrix *g;
    bool conv;

public:
    sta();
//    sta(const sta& orig);
    virtual ~sta();
    double time() const;
    void op();
    void print() const;
private:
    void init();
    void load();
    int iterate(const int& , const int&);
    void solve();
    bool convergence();
    void printSolution();
};

#endif /* STA_H */

