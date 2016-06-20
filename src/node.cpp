#include "global.h"
#include "equation.h"
#include "model.h"
#include "matrix.h"
#include "node.h"

node::node() : node("gnd") {}

node::node(string s) {
    name = s;
    value[0] = 0.0;
    value[1] = 0.0;
    valueDelta = 0.0;
    icFlag = false;
    icValue = 0.0;
    trValue[0] = 0.0;
    trValue[1] = 0.0;
    trValue[2] = 0.0;
    trValue[3] = 0.0;
    dotValue = 0.0;
    trDotValue[0] = 0.0;
    trDotValue[1] = 0.0;
    ncon = 0;
    if (name == "gnd" || name == "GND" || name == "0")
        ground = true;
    else
        ground = false;
    printFlag = false;
    model::currentModel->add(this);
}

node::node(string s, double d) : node(s) {
    icValue = d;
}

//node::node(node& n) {
//    name = n.name;
//    value = 0.0;
//    trValue[0] = 0.0;
//    trValue[1] = 0.0;
//    trValue[2] = 0.0;
//    trValue[3] = 0.0;
//    dotValue = 0.0;
//    trDotValue[0] = 0.0;
//    trDotValue[1] = 0.0;
//    ground = n.ground;
//    printFlag = n.printFlag;
//}

ostream& operator <<(ostream& os, const node& n) {
    os << "node " << &n << " = {";
    os << " name=\"" << n.name << "\"";
    os << ", value[0]=" << n.value[0];
    os << ", index=" << n.index;
    os << " }";
    return(os);
}

void node::reset() {
    value[0] = 0.0;
    value[1] = 0.0;
    dotValue = 0.0;
    trDotValue[0] = 0.0;
    trDotValue[1] = 0.0;
}

void node::init() {
}

void node::addConnection() {
    ++ncon;
}
