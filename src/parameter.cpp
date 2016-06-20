#include "parameter.h"
#include "model.h"
#include "global.h"

parameter_::parameter_() {
    name = "";
    model::currentModel->add(this);
}

parameter_::parameter_(const string s) {
    name = s;
    model::currentModel->add(this);
}

ostream& operator <<(ostream& os, const parameter_& p) {
    return os;
}


