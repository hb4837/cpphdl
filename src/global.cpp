#include "global.h"

const double inf = HUGE_VAL;
const double undef = HUGE_VAL - 1;

double gmin = 1.0e-12;
double temp = 27.0;
double tempk = 300.15;

double time() {
    return sim::getInstance()->time();
}

