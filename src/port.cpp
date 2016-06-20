#include "port.h"
#include "model.h"

port::port(string s) {
    name = s;
    nodePtr = 0;
    model::add(this);
}

port::port() {
    nodePtr = 0;
    model::add(this);
}

port::port(node& n) {
    nodePtr = &n;
    model::add(this);
}

void port::operator() (node& n) {
    nodePtr = &n;
}

ostream& operator <<(ostream& os, const port& p) {
    os << "port " << &p << " = {";
    os << " name=\"" << p.name << "\"";
    os << " nodePtr=" << p.nodePtr;
    os << " }";
    return os;
}