#include "indent_facet.hpp"
#include "model.h"
#include "node.h"
#include "quantity.h"
#include "equation.h"

list<model *> model::instances = (list<model *>) 0;
model *topLevel = new model("topLevel");
model *model::currentModel = topLevel;

model::model(const string s) {
    name = s;
    model::instances.push_back(this);
    currentModel = this;
}

model::model(const char *s) {
    name = string(s);
    model::instances.push_back(this);
    currentModel = this;
}

void model::add(node *n) {
    currentModel->nodes.push_back(n);
}

void model::add(quantity_ *q) {
    q->modelPtr = currentModel;
    q->name += "(" + q->modelPtr->name + ")";
    currentModel->quantities.push_back(q);
}

void model::add(equation *e) {
    currentModel->equationsList.push_back(e);
    e->modelPtr = currentModel;
    e->name += "f" + to_string(currentModel->equationsList.size()) + "(" + e->modelPtr->name + ")";
}

void model::add(parameter_ *p) {
    currentModel->params.push_back(p);
}

void model::add(port *p) {
    currentModel->ports.push_back(p);
}

void model::init() {
    //    for(auto ieq = equationsList.begin(); ieq != equationsList.end(); ++ieq) {
    //        (*ieq)->init();
    //    }

    //    for(auto iq = quantities.begin(); iq != quantities.end(); ++iq) {
    //        quantity_ *q = (*iq);
    //        if (q->isAcrossQuantity() || q->isThroughQuantity()) {
    //            cerr << "name=" << q->name << endl;
    //            cerr << "from=" << q->getFromPort() << " -> " << q->getFromNode() << endl;
    //            cerr << "to=" << q->getToPort() << " -> " << q->getToNode() << endl;
    //        }
    //    }
}

//void model::load() {
//    list<equation *>::iterator ieq;
//
//    for(ieq = equations.begin(); ieq != equations.end(); ++ieq) {
//        (*ieq)->load();
//    }
//}

void model::parameters() {
}

//void model::evaluate() {
//#ifdef DEBUG
//    cout << endl << endl;
//    cout << "------------------------" << endl;
//    cout << name << "->evaluate()" << endl;
//    cout << "------------------------" << endl;
//    cout << endl;
//#endif
//
////    curEquation = equations.begin();
////    behaviour();
//}

void model::equations() {
}

void model::limit() {
}

void model::load() {
}

ostream& operator<<(ostream& os, const model& m) {
    os << "model " << &m << " = {" << endl;
    os << im::push;
    os << "name=\"" << m.name << "\"" << endl;

    os << "ports={";
    for (auto ip = m.ports.begin(); ip != m.ports.end(); ++ip) {
        if (ip != m.ports.begin())
            os << ", ";
        os << *(*ip);
    }
    os << "}" << endl;

    os << "nodes={";
    for (auto in = m.nodes.begin(); in != m.nodes.end(); ++in) {
        if (in != m.nodes.begin())
            os << ", ";
        os << *(*in);
    }
    os << "}" << endl;

    os << "equations={";
    for (auto ieq = m.equationsList.begin(); ieq != m.equationsList.end(); ++ieq) {
        if (ieq != m.equationsList.begin())
            os << endl;
        os << *ieq;
    }
    os << "}" << endl;

    os << "quantities = {";
    for (auto iq = m.quantities.begin(); iq != m.quantities.end(); ++iq) {
        if (iq != m.quantities.begin())
            os << ", ";
        os << *(*iq);
    }
    os << "}" << endl;
    cout << im::pop;
    os << "}" << endl;
    return (os);
}

//void model::getSolution() {
//}

void model::printInstances() {

    for (auto iq = model::instances.begin(); iq != model::instances.end(); ++iq) {
        model *m = *iq;
        cout << "model " << m << " {" << im::push << endl;
        cout << "name=\"" << m->name << "\"" << endl;
        cout << "nodes = {" << im::push << endl;
        for (auto in = m->nodes.begin(); in != m->nodes.end(); ++in)
            cout << *(*in) << endl;
        cout << im::pop << "}" << endl;
        cout << "ports = {" << im::push << endl;
        for (auto ip = m->ports.begin(); ip != m->ports.end(); ++ip)
            cout << *(*ip) << endl;
        cout << im::pop << "}" << endl;
        cout << "nodes = {" << im::push << endl;
        for (auto ip = m->nodes.begin(); ip != m->nodes.end(); ++ip)
            cout << *(*ip) << endl;
        cout << im::pop << "}" << endl;
        cout << "quantities = {" << im::push << endl;
        for (auto ip = m->quantities.begin(); ip != m->quantities.end(); ++ip)
            cout << *(*ip) << endl;
        cout << im::pop << "}" << endl;
        cout << im::pop << "}" << endl;
    }
}
