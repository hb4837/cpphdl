#include "spice.h"
#include "pll.h"

int main(int argc, char *argv[]) {
    NODE(gnd);
    NODE(in);
    NODE(out);

    vco vco1("vco1", in, gnd, out, gnd, 1.0, 7.0e6, 30.0e3);
    r r1 ("r1", out, gnd, 1.0e6);
    vdc vin("vin", in, gnd, 1.0);

    op::monitor = 1;
    kernel::op();

//  plot::output = *argv;
//  plot::add("in", "pll_in", "lf_out", "vco_out");
//
//  tran::monitor = 1;
//  tran::tsmin = 1.0e-21;
//  tran(1.0e-7, 300.0e-6, 0.0, 1.0e-9);
}
