#include "spice.h"

int main(int argc, char *argv[]) {    
    NODE(gnd);
    NODE(n1);
    NODE(n2);

    v v1("vin");
    v1.pos(n1);
    v1.neg(gnd);
    v1.dc(3.3);
    // v1.pulse(0.0, 1.0, 1.0e-6, 1.0e-9, 1.0e-9, 10.0e-6);

    r r1("r1");
    r1.pos(n1);
    r1.neg(n2);
    r1.rnom(1e3);

    r r2("r2");
    r2.pos(n2);
    r2.neg(gnd);
    r2.rnom(1e3);

//    c c1("c1");
//    c1.pos(n2);
//    c1.neg(gnd);
//    c1.value(1e-9);

    sta solver;
    op::monitor=1;
    
    // model::printInstances();
    
//    op::monitor = 1;
    solver.op();
    solver.print();

//  tran::monitor = 10;
//  tran(1.0e-10, 21e-6, 0.0, 1.0e-7);
//
//  plot::output = *argv;
//  plot::add("n1", "n2");
}
