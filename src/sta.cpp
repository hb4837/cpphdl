#include "sta.h"
#include "util.h"
#include "indent_facet.hpp"

int op::monitor = 0;
int op::newton::itlim = 100;
double op::newton::abstol = 1.0e-12;
double op::newton::reltol = 1.0e-3;
double op::newton::vntol = 1.0e-6;

sta::sta() {
    t = 0.0;
    g = 0;
    numUnknowns = 0;
    numNodes = 0;
    sim::setInstance(this);
    conv = false;
}

//sta::sta(const sta& orig) {
//}

sta::~sta() {
}

double sta::time() const {
    return t;
}

void sta::init() {

    // Collect all nodes
#ifdef DEBUG
    cout << "Nodes:" << endl;
#endif
    for (auto im = model::getInstances().begin(); im != model::getInstances().end(); ++im) {
        model *m = *im;

        m->init();

        const list<port *> ports = m->getPorts();
        for (auto ip = m->getPorts().begin(); ip != m->getPorts().end(); ++ip) {
            node *n = (*ip)->getNode();
            if (!n->isGround()) {
                string key = n->getName();
                if (nodes.find(key) == nodes.end()) {
                    nodes[key] = n;
#ifdef DEBUG
                    cout << numUnknowns << ": " << n->getName() << endl;
#endif
                    n->setIndex(numUnknowns++);
                    numNodes++;
                }
            }
        }
    }

    int i = numNodes;
#ifdef DEBUG
    cout << "Equations:" << endl;
#endif
    for (auto im = model::getInstances().begin(); im != model::getInstances().end(); ++im) {
        model *m = *im;
        const list<equation *>& equations = m->getEquations();
        for (auto ieq = equations.begin(); ieq != equations.end(); ++ieq) {
#ifdef DEBUG
            cout << i << ": " << (*ieq)->getName() << endl;
#endif
            (*ieq)->setIndex(i++);
        }
    }

#ifdef DEBUG
    cout << "Quantities:" << endl;
#endif
    // Collect all through quantities
    for (auto im = model::getInstances().begin(); im != model::getInstances().end(); ++im) {
        model *m = *im;

        for (auto iq = m->getQuantities().begin(); iq != m->getQuantities().end(); ++iq) {
            quantity_ *q = *iq;
            if (q->isThroughQuantity()) {
#ifdef DEBUG
                cout << numUnknowns << ": " << q->getName() << endl;
#endif
                q->setIndex(numUnknowns++);

                string key = q->getName();
                if (currents.find(key) == currents.end()) {
                    currents[key] = q;
                }
            }
        }
    }

    g = new matrix(numUnknowns);
}

void sta::load() {

    g->clear();

    for (auto itm = model::getInstances().begin(); itm != model::getInstances().end(); ++itm) {

        // Execute model equations
        model *m = *itm;
        m->load();
        const list<equation *>& equations = m->getEquations();

        // Incidence matrix A
        for (auto ieq = equations.begin(); ieq != equations.end(); ++ieq) {
            equation *eq = *ieq;

            for (auto iq = eq->getQuantities().begin(); iq != eq->getQuantities().end(); ++iq) {
                quantity_ *q = (*iq);

                node *from = q->getFromNode();
                if (!from->isGround()) {
                    if (!q->isAcrossQuantity())
                        g->addEntry(from->getIndex(), q->getIndex(), 1.0);
                }

                node *to = q->getToNode();
                if (!to->isGround()) {
                    if (!q->isAcrossQuantity())
                        g->addEntry(to->getIndex(), q->getIndex(), -1.0);
                }
            }
        }

        // df/dv and df/di
        for (auto ieq = equations.begin(); ieq != equations.end(); ++ieq) {
            equation *eq = *ieq;

            // df/dv
            for (auto in = nodes.begin(); in != nodes.end(); ++in) {
                node *n = (*in).second;

                for (auto iq = eq->getQuantities().begin(); iq != eq->getQuantities().end(); ++iq) {
                    quantity_ *q = *iq;

                    if (q->isAcrossQuantity()) {
                        if (q->getFromNode() == n)
                            g->addEntry(eq->getIndex(), n->getIndex(), eq->getJ(*q));
                        if (q->getToNode() == n)
                            g->addEntry(eq->getIndex(), n->getIndex(), -eq->getJ(*q));
                    }
                }
            }

            // df/di
            for (auto iq = eq->getQuantities().begin(); iq != eq->getQuantities().end(); ++iq) {
                quantity_ *q = *iq;

                if (q->isThroughQuantity()) {
                    g->addEntry(eq->getIndex(), q->getIndex(), eq->getJ(*q));
                }
            }

            // -f
            g->addRhsEntry(eq->getIndex(), -eq->getF());
        }

    }
}

void sta::op() {
    init();

    cout << "---------------" << endl;
    cout << "Operating point" << endl;
    cout << "---------------" << endl;
    cout << endl;

    int iterno = iterate(op::newton::itlim, op::monitor);
    //  if(icFlag) {
    //      icFlag = false;
    //      iterno = iterate(op::newton::itlim, op::monitor);
    //  }

    //  if(!conv && op::ramp::src::itlim != 0) {
    //      cout << "*op*: newton algorithm failed to converge" << endl;
    //      cout << "      trying source stepping..." << endl;
    //      iterno = sourceStepping(op::ramp::src::itlim);
    //      if(conv && iterno <= op::ramp::src::itlim)
    //          cout << "*message*: source stepping algorithm succeeded." << endl;
    //      else
    //          error("source stepping algorithm failed.");
    //  }

    if (!conv) {
        cout << "*op*: no convergence in op analysis" << endl;
        cout << "      last solution vector:" << endl;
        // ...
        // Write solution to nodes & quantities

        printSolution();
        exit(1);
    }

    //  if(analysis == OP) {
    //      cout << endl << "*op*: operating point analysis results:" << endl;
    //      for(in = nodes.begin(); in != nodes.end(); ++in) {
    //          if((*in)->ground)
    //              continue;
    //          cout << setw(15) << ((*in)->name).data();
    //          cout << setw(15) << (*in)->v0 << endl;
    //      }
    //      //for(iq = quantities.begin(); iq != quantities.end(); ++iq) {
    //          //cout << setw(15) << ((*iq)->name).data();
    //          //cout << setw(15) << (*iq)->v0 << endl;
    //      //}
    //  }
    //  status = OP_DONE;

    cout << endl;
    cout << "Operating point analysis finished" << endl;
    cout << endl;
    cout << "Statistics: " << endl;
    cout << "    number of iterations:   " << setw(15) << iterno << endl;
    cout << "    number of unknowns:     " << setw(15) << g->getSize() << endl;
    cout << endl;

}

int sta::iterate(const int& itlim, const int& mon) {
    int iterno = 0;

    if (mon > 0) {
        cout << endl;
        cout << setw(5) << "it";
        cout << setw(5) << "+/-";
        cout << endl;
    }

    do {
#ifdef DEBUG
        cout << endl << endl << endl;
        cout << "===========================" << endl;
        cout << "Starting iteration no." << setw(5) << iterno + 1 << endl;
        cout << "===========================" << endl;
#endif
        load();
        solve();

        conv = convergence();

        iterno++;

        if (mon > 0) {
            if (iterno % mon == 0) {
                cout << setw(5) << iterno;
                if (conv)
                    cout << setw(5) << "+";
                else {
                    //                    if(limited) {
                    //                        cout << setw(5) << "*";
                    //                        limited = false;
                    //                    }
                    //                    else {
                    cout << setw(5) << "-";
                    //                    }
                }
                cout << endl;
            }
        }
        //        else
        //            limited = false;
    } while (!conv && iterno <= itlim);

    return (iterno);
}

void sta::solve() {
    // g->info();
    g->factor();
    g->solve();
    // g->info();

    // Write solution to nodes & quantities
    cout << "[";
    for (auto in = nodes.begin(); in != nodes.end(); ++in) {
        node *n = (*in).second;
        double delta = g->getSolution(n->getIndex());
        n->pushAddValue(delta);
        if (in != nodes.begin())
            cout << ", ";
        cout << n->getValue();
    }
    cout << "]";

    cout << "[";
    for (auto ic = currents.begin(); ic != currents.end(); ++ic) {
        quantity_ *q = (*ic).second;

        // if (q->isThroughQuantity()) {
        double delta = g->getSolution(q->getIndex());
        q->pushAddValue(delta);
        if (ic != currents.begin())
            cout << ", ";
        cout << q->getValue();
        // }
    }
    cout << "]" << endl;
}

void sta::print() const {
    for (auto in = nodes.begin(); in != nodes.end(); ++in) {
        cout << (*in).second->getName() << ": " << (*in).second->getValue() << endl;
    }
    for (auto ic = currents.begin(); ic != currents.end(); ++ic) {
        cout << (*ic).second->getName() << ": " << (*ic).second->getValue() << endl;
    }
}

void sta::printSolution() {
    cout << "[";
    for (auto in = nodes.begin(); in != nodes.end(); ++in) {
        node *n = (*in).second;
        double delta = g->getSolution(n->getIndex());
        n->pushAddValue(delta);
        if (in != nodes.begin())
            cout << ", ";
        cout << n->getValue();
    }
    cout << "] [";
    for (auto ic = currents.begin(); ic != currents.end(); ++ic) {
        quantity_ *q = (*ic).second;

        double delta = g->getSolution(q->getIndex());
        q->pushAddValue(delta);
        if (ic != currents.begin())
            cout << ", ";
        cout << q->getValue();
    }
    cout << "]" << endl;
}

bool sta::convergence() {
    for (auto in = nodes.begin(); in != nodes.end(); ++in) {
        node *n = (*in).second;
        double v0 = n->getValue(0);
        double v1 = n->getValue(1);
        double absVal0 = fabs(v0);
        double absVal1 = fabs(v1);
        double errVal = fabs(v0 - v1);
        double absMin = (absVal0 < absVal1) ? absVal0 : absVal1;
        double limit = op::newton::vntol + op::newton::reltol * absMin;
        if (errVal > limit)
            return false;
    }

    for (auto ic = currents.begin(); ic != currents.end(); ++ic) {
        quantity_ *q = (*ic).second;
        double v0 = q->getValue(0);
        double v1 = q->getValue(1);
        double absVal0 = fabs(v0);
        double absVal1 = fabs(v1);
        double errVal = fabs(v0 - v1);
        double absMin = (absVal0 < absVal1) ? absVal0 : absVal1;
        double limit = op::newton::abstol + op::newton::reltol * absMin;
        if (errVal > limit)
            return false;
    }

    return (true);
}
