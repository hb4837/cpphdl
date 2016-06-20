#include "spice.h"

int main(int argc, char *argv[]) {
    node gnd;
    node n1("n1");
    node n2("n2");

    vdc v1("v1", n1, gnd, 1000.0);
    d d1("d1", n1, n2);
    //d1.bv = 5.0;
    r r1("r1", n2, gnd, 1.0e6);

//    //op::newton::itlim = 3;
//    op::newton::reltol = 1.0e-9;
//    op::monitor = 1;
//    op();

    // plot::output = *argv;
    // plot::add("vd(d1)", "id(d1)");
    // kernel::dc(v1, -5.7, 0.7, 0.01);
}
