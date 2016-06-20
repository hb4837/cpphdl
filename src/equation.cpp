#include <sstream>
#include "equation.h"
#include "matrix.h"
#include "node.h"
#include "quantity.h"
#include "model.h"
#include "indent_facet.hpp"

equation::equation() {
    name = "";
    type = NONE;
    index = 0;
    f = 0.0;
    model::add(this);
}

equation::equation(quantity_& q1) : equation() {
    quantities.push_back(&q1);
}

equation::equation(quantity_& q1, quantity_& q2) : equation() {
    quantities.push_back(&q1);
    quantities.push_back(&q2);
}

equation::~equation() {
    //cerr << "equation::~equation(): name=" << name << endl;
}

void equation::setup() {
    //    list<quantity_ *>::const_iterator iq;
    //
    //    for(iq = modelPtr->quantities.begin(); iq != modelPtr->quantities.end(); ++iq) 
    //        add(*iq);
}

ostream& operator<<(ostream& os, const equation& e) {
    os << "equation " << &e << " = {";
    os << im::push;
    os << endl;
    os << "name = \"" << e.name << "\"" << endl;
    os << "index = " << e.index << endl;
    os << "f = " << e.f << endl;
    os << "quantities = {";
    os << im::push;
    for (auto iq = e.quantities.begin(); iq != e.quantities.end(); ++iq) {
        if (iq != e.quantities.begin())
            os << ", ";
        os << (*iq)->getName();
    }
    cout << im::pop;
    os << "}" << endl;
    os << "j = {";
    for (auto ij = e.j.begin(); ij != e.j.end(); ++ij) {
        if (ij != e.j.begin())
            os << ", ";
        os << " \"" << (*ij).first->getName() << "\"->" << (*ij).second;
    }
    os << " }" << endl;
    cout << im::pop;
    os << "}";
    return (os);
}

//void equation::operator =(const quantity_& q) {
//    f = q.value[0];
//    quantities = q.quantities;
//    deriv = q.grad;
//}

void equation::init() {
    //list<quantity_ *>::const_iterator iq1;
    //map<const string, quantity_ *>::const_iterator iq2;

    //model *m = kernel::getCurModel();

    //for(iq1 = m->quantities.begin(); iq1 != m->quantities.end(); ++iq1) {
    //for(iq2 = quantities.begin(); iq2 != quantities.end(); ++iq2) {
    //}
    //}
}

void equation::load() {
}

