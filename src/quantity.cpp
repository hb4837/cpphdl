#include "quantity.h"
#include "node.h"
#include "equation.h"
#include "matrix.h"
#include "port.h"
#include "indent_facet.hpp"

#undef EQNDEBUG

quantity_::quantity_() : quantity_("") {
}

quantity_::quantity_(const string& s) {
    name = s;
    type = UNKNOWN;
    index = 0;
    value[0] = 0.0;
    value[1] = 0.0;
    dotValue = 0.0;
    trDotValue[0] = 0.0;
    trDotValue[1] = 0.0;
    icFlag = false;
    icValue = 0.0;
    valueDelta = 0.0;
    limitFlag = false;
    trValue[0] = 0.0;
    trValue[1] = 0.0;
    trValue[2] = 0.0;
    trValue[3] = 0.0;
    printFlag = false;
    dynamicFlag = false;
    implicitFlag = true;
    modelPtr = 0;
    equationPtr = 0;
}

quantity_::~quantity_() {
}

ostream& operator<<(ostream& os, const quantity_& q) {
    string typeName;
        
    switch (q.type) {
        case(quantity_::UNKNOWN): typeName = "UNKNOWN";
            break;
        case(quantity_::ACROSS): typeName = "ACROSS";
            break;
        case(quantity_::THROUGH): typeName = "THROUGH";
            break;
        case(quantity_::DOUBLE): typeName = "DOUBLE";
            break;
        case(quantity_::INTERMEDIATE): typeName = "INTERMEDIATE";
            break;
    }

    os << "quantity_ " << &q << " = {" ;
    os << im::push;
    os << endl;
    
    os << "name=\"" << q.name << "\"" << endl;
    // os << " type=" << typeName;
    //    os << "index = " << q.index;
    //    os << "value = " << q.value;
    //    os << "quantities = {";
    //    for(iq = q.quantities.begin(); iq != q.quantities.end(); ++iq)
    //        os << " " << iq->second->name;
    //    os << "}";
    //    for(iq = q.quantities.begin(); iq != q.quantities.end(); ++iq) {
    //        os << "grad[" << iq->second->getName() << "] = "
    //           << iq->second->getGrad(iq->second->getName()) << endl;
    //    }
    //    os << "trValue = {" << q.trValue[0] << ", " << 
    //        q.trValue[1] << ", " << q.trValue[2] << ", " << 
    //        q.trValue[3] << "}" << endl;
    //    os << "value[1] = " << q.value[1] << endl;
    //    os << "delta = " << q.valueDelta << endl;
    //    os << "dotValue = " << q.dotValue << endl;
    //    os << "trDotValue = {" << q.trDotValue[0] << ", " << 
    //        q.trDotValue[1] << "}" << endl;
    //    os << "printFlag = " << 
    //        (q.printFlag  == true ? "true" : "false") << endl;
    //    os << "dynamicFlag = " << 
    //        (q.dynamicFlag  == true ? "true" : "false") << endl;
    //    os << "implicitFlag = " << 
    //        (q.implicitFlag  == true ? "true" : "false") << endl;
    //    //os << indent << "quantities = ";
    //    //for(iq = q.quantities.begin(); iq != q.quantities.end(); ++iq) {
    //        //os << " " << ((*iq)->name).data();
    //    //}
    //    //os << endl;
    cout << im::pop;
    os << "}";
    return (os);
}

void quantity_::reset() {
    value[0] = 0.0;
    value[1] = 0.0;
    dotValue = 0.0;
    trDotValue[0] = 0.0;
    trDotValue[1] = 0.0;
}

double double_cast(const quantity_& q) {
    return q.value[0];
}

//template<>
//quantity<across>::quantity() : quantity_() {
//    type = ACROSS;
//    from = 0;
//    to = 0;
//    implicitFlag = false;
//    model::add(this);
//}
//
//template<>
//quantity<across>::quantity(const string& s) : quantity_(s) {
//    type = ACROSS;
//    from = 0;
//    to = 0;
//    implicitFlag = false;
//    model::add(this);
//}

template<>
void quantity<across>::startvalue(const double& d) {
    icFlag = true;
    icValue = d;
}

//template<>
//quantity<through>::quantity() : quantity_() {
//    type = THROUGH;
//    from = 0;
//    to = 0;
//    model::add(this);
//}
//
//template<>
//quantity<through>::quantity(const string& s) : quantity_(s) {
//    type = THROUGH;
//    from = 0;
//    to = 0;
//    model::add(this);
//}

template<>
quantity<double>::quantity() : quantity_() {
    type = DOUBLE;
    from = 0;
    to = 0;
    implicitFlag = false;
    model::add(this);
}

template<>
quantity<double>::quantity(const string& s) : quantity_(s) {
    type = DOUBLE;
    from = 0;
    to = 0;
    quantities.push_back(this);
    grad.push_back(1.0);
    model::add(this);
}

template<>
quantity<across>::quantity(const string& s, port& f, port& t)
: quantity_(s), from(&f), to(&t) {
    type = ACROSS;
    quantities.push_back(this);
    grad.push_back(1.0);
    from = &f;
    to = &t;
    implicitFlag = false;
    model::add(this);
}

template<>
quantity<through>::quantity(const string& s, port& f, port& t)
: quantity_(s), from(&f), to(&t) {
    type = THROUGH;
    quantities.push_back(this);
    grad.push_back(1.0);
    from = &f;
    to = &t;
    model::add(this);
}

template<>
quantity<through>::quantity(const string& s, port& f)
: quantity_(s), from(&f) {
    type = THROUGH;
    model::add(this);
    from = &f;
}

template<>
void quantity<across>::operator=(const double& d) {
    value[0] = d;
}

//void quantity<double>::operator =(const quantity_& q) {
//value = q.value;
//}

template<>
quantity<double>::quantity(quantity_ q) : quantity_(q) {
    //cerr << "quantity(" << q << ")" << endl;
    //value[0] = q.getValue();
    //value[1] = q.value[1];
    //dotValue = q.dotValue;
    //trDotValue[0] = q.trDotValue[0];
    //trDotValue[1] = q.trDotValue[1];
}

quantity_ quantity_::dot() {
    quantity_ r;
    return r;
}

//quantity_ pow(const quantity_& q, const double& d) {
//#ifdef EQNDEBUG
//    string f = "pow(" + q.name + ", (double))";
//    quantity_ r(f);
//#else
//    quantity_ r;
//#endif
//    map<const string, quantity_ *>::iterator iq;
//    double dq_dx;
//
//    r.type = quantity_::INTERMEDIATE;
//
//    r.add(q);
//
//    r.value[0] = pow(q.value[0], d);
//
//    double df_dq = d * pow(q.value[0], d - 1.0);
//    for (iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//        dq_dx = q.getGrad(iq->second->name);
//        r.grad[iq->second->name] = df_dq * dq_dx;
//    }
//#ifdef EQNDEBUG
//    cout << r << endl;
//    r.check();
//#endif
//    return (r);
//}

//quantity_ exp(const quantity_& q) {
//#ifdef EQNDEBUG
//    string f = "exp(" + q.name + ")";
//    quantity_ r(f);
//#else
//    quantity_ r;
//#endif
//    map<const string, quantity_ *>::iterator iq;
//    double g;
//
//    r.type = quantity_::INTERMEDIATE;
//
//    r.add(q);
//
//    r.value[0] = exp(q.value[0]);
//    for (iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//        g = q.getGrad(iq->second->name);
//        r.grad[iq->second->name] = g * r.value[0];
//    }
//#ifdef EQNDEBUG
//    cout << r << endl;
//    r.check();
//#endif
//    return (r);
//}

//quantity_ limexp(const quantity_& q) {
//#ifdef EQNDEBUG
//    string f = "limexp(" + q.name + ")";
//    quantity_ r(f);
//#else
//    quantity_ r;
//#endif
//    map<const string, quantity_ *>::iterator iq;
//    double g;
//    double arg;
//
//    r.type = quantity_::INTERMEDIATE;
//
//    r.add(q);
//
//#define EXP_MAXARG (50.0)
//    arg = exp(EXP_MAXARG);
//
//    if (q.value[0] > EXP_MAXARG) {
//        r.value[0] = arg * (1.0 + q.value[0] - EXP_MAXARG);
//        for (iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//            g = q.getGrad(iq->second->name);
//            r.grad[iq->second->name] = g * arg;
//        }
//    } else {
//        r.value[0] = exp(q.value[0]);
//        for (iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//            g = q.getGrad(iq->second->name);
//            r.grad[iq->second->name] = g * r.value[0];
//        }
//    }
//#ifdef EQNDEBUG
//    cout << r << endl;
//    r.check();
//#endif
//    return (r);
//}

//quantity_ log(const quantity_& q) {
//#ifdef EQNDEBUG
//    string f = "log(" + q.name + ")";
//    quantity_ r(f);
//#else
//    quantity_ r;
//#endif
//    map<const string, quantity_ *>::iterator iq;
//    double g;
//
//    r.type = quantity_::INTERMEDIATE;
//
//    r.add(q);
//    r.value[0] = log(q.value[0]);
//    for (iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//        g = q.getGrad(iq->second->name);
//        r.grad[iq->second->name] = g / q.value[0];
//    }
//#ifdef EQNDEBUG
//    cout << r << endl;
//    r.check();
//#endif
//    return (r);
//}

quantity_ operator-(const quantity_& q) {
#ifdef EQNDEBUG
    string f = "-" + q.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = -q.value[0];
    
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        double g = q.grad[i];
        r.grad[i] = -g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator+(const quantity_& q1, const quantity_& q2) {
#ifdef EQNDEBUG
    string f = q1.name + "+" + q2.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g1, g2;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q1);
    r.quantities.push_back((quantity_ *) &q2);

    r.value[0] = q1.value[0] + q2.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g1 = q1.grad[i];
        g2 = q2.grad[i];
        r.grad[i] = g1 + g2;
    }

#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator+(const quantity_& q, const double& d) {
#ifdef EQNDEBUG
    string f = q.name + "+ (double)";
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = q.value[0] + d;
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator+(const double& d, const quantity_& q) {
#ifdef EQNDEBUG
    string f = "(double)+" + q.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = d + q.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return r;
}

quantity_ operator-(const quantity_& q1, const quantity_& q2) {
#ifdef EQNDEBUG
    string f = q1.name + "-" + q2.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g1, g2;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q1);
    r.quantities.push_back((quantity_ *) &q2);

    r.value[0] = q1.value[0] - q2.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g1 = q1.grad[i];
        g2 = q2.grad[i];
        r.grad[i] = g1 - g2;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator-(const quantity_& q, const double& d) {
#ifdef EQNDEBUG
    string f = q.name + "-(double)";
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);
    
    r.value[0] = q.value[0] - d;
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator-(const double& d, const quantity_& q) {
#ifdef EQNDEBUG
    string f = "(double)-" + q.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g;
    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = d - q.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = -g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator*(const quantity_& q1, const quantity_& q2) {
#ifdef EQNDEBUG
    string f = q1.name + "*" + q2.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g1, g2;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q1);
    r.quantities.push_back((quantity_ *) &q2);

    r.value[0] = q1.value[0] * q2.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g1 = q1.grad[i];
        g2 = q2.grad[i];
        r.grad[i] = g1 * q2.value[0] + q1.value[0] * g2;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator*(const quantity_& q, const double& d) {
#ifdef EQNDEBUG
    string f = q.name + "*(double)";
    quantity_ r(f);
#else
    quantity_ r;
#endif
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = q.value[0] * d;
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = g * d;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator*(const double& d, const quantity_& q) {
#ifdef EQNDEBUG
    string f = "(double)*" + q.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    map<const string, quantity_ *>::iterator iq;
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = d * q.value[0];

    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = g * d;
    }

    //???
    //if(time_domain()) {
    //r.trValue[0] = d * q.trValue[0];
    //r.trValue[1] = d * q.trValue[1];
    //r.dotValue = d * q.dotValue;
    //r.trDotValue[0] = d * q.trDotValue[0];
    //r.trDotValue[1] = d * q.trDotValue[1];
    //cerr << "!!" << r.getName() << ": " << r.trDotValue[0] << endl;
    //}

#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator/(const quantity_& q1, const quantity_& q2) {
#ifdef EQNDEBUG
    string f = q1.name + "/" + q2.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    map<const string, quantity_ *>::iterator iq;
    double g1, g2;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q1);
    r.quantities.push_back((quantity_ *) &q2);

    r.value[0] = q1.value[0] / q2.value[0];
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g1 = q1.grad[i];
        g2 = q2.grad[i];
        r.grad[i] = (g1 * q2.value[0] - g2 * q1.value[0]) / (q2.value[0] * q2.value[0]);
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator/(const quantity_& q, const double& d) {
#ifdef EQNDEBUG
    string f = q.name + "/(double)";
    quantity_ r(f);
#else
    quantity_ r;
#endif
    map<const string, quantity_ *>::iterator iq;
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = q.value[0] / d;
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad.push_back(g / d);
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

quantity_ operator/(const double& d, const quantity_& q) {
#ifdef EQNDEBUG
    string f = "(double)/" + q.name;
    quantity_ r(f);
#else
    quantity_ r;
#endif
    map<const string, quantity_ *>::iterator iq;
    double g;

    r.type = quantity_::INTERMEDIATE;

    r.quantities.push_back((quantity_ *) &q);

    r.value[0] = d / q.value[0];
    double tmp = -d / (q.value[0] * q.value[0]); // ???
    for (unsigned int i = 0; i < r.quantities.size(); ++i) {
        g = q.grad[i];
        r.grad[i] = tmp * g;
    }
#ifdef EQNDEBUG
    cout << r << endl;
    r.check();
#endif
    return (r);
}

void operator==(const quantity_& q1, const quantity_& q2) {
//#ifdef EQNDEBUG
//    string fname = q1.name + "==" + q2.name;
//    quantity_ f(fname);
//#else
//    quantity_ f;
//#endif
//
//    equation *eq;
//    if (q1.getEquation() == 0) {
//        eq = new equation();
//        model::add(eq);
//    } else {
//        eq = q1.getEquation();
//    }
//
//    f.value[0] = q1.value[0] - q2.value[0];
//    
//    for (unsigned int i = 0; i < q1.quantities.size(); ++i) {
//        f.quantities.push_back(q1.quantities[i]);
//    }
//    
//    for (unsigned int i = 0; i < q2.quantities.size(); ++i) {
//        quantity_ *qi = q2.quantities[i];
//        
//        bool found = false;
//        for (unsigned int j = 0; j < f.quantities.size(); ++j) {
//            quantity_ *qj = f.quantities[j];
//            if (qi == qj)
//                found = true;
//        }
//        if (!found)
//            f.quantities.push_back(q1.quantities[i]);
//    }
//    
//    for (unsigned int i = 0; i < f.quantities.size(); ++i) {
//        f.grad.push_back(q1.grad[i] - q2.grad[i]);
//    }
//    
//    eq->f = f.value[0];
//    eq->quantities = f.quantities;
//    eq->j = f.j;
}

void operator==(const quantity_& q1, const double& d) {
//#ifdef EQNDEBUG
//    string fname = q1.name + "==" + to_string(d);
//    quantity_ f(fname);
//#else
//    quantity_ f;
//#endif
//
//    // Create equation object if it does not exist yet
//    equation *eq;
//    if (q1.getEquation() == 0) {
//        eq = new equation();
//        model::add(eq);
//    } else {
//        eq = q1.getEquation();
//    }
//
//    // Store result value
//    f.value[0] = q1.value[0] - d;
//    for (unsigned int i = 0; i < f.quantities.size(); ++i) {
//        f.grad[i] = q1.grad[i];
//    }
//    
//    eq->f = f.value[0];
//    eq->quantities = f.quantities;
//    eq->j = f.grad;
//    
//    //#ifdef EQNDEBUG
//    //    string f = q.name + "==(double)";
//    //    quantity_ r(f);
//    //#else
//    //    quantity_ r;
//    //#endif
//    //    map<const string, quantity_ *>::const_iterator iq;
//    //    double g;
//    //    equation *eq;
//    //    
//    //    r.type = quantity_::INTERMEDIATE;
//    //
//    //    r.add(q);
//    //
//    //    if(kernel::initStatus()) {
//    //        eq = new equation();
//    //        kernel::getCurModel()->addEquation(eq);
//    //
//    //        for(iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//    //            if(iq->second->type == quantity_::THROUGH)
//    //               //iq->second->type == quantity_::DOUBLE)
//    //            {
//    //                iq->second->setExplicit();
//    //            }
//    //        }
//    //    }
//    //    else {
//    //        eq = kernel::getCurModel()->getEquation();
//    //
//    //        r.value[0] = q.value[0] - d;
//    //        for(iq = r.quantities.begin(); iq != r.quantities.end(); ++iq) {
//    //            g = q.getGrad(iq->second->name);
//    //            r.grad[iq->second->name] = g;
//    //        }
//    //        *eq = r;
//    //    }
}

void quantity_::operator=(const quantity_& q) {
    setExplicit();
    value[1] = value[0];
    value[0] = q.value[0];
    dotValue = q.dotValue;
    valueDelta = value[0] - value[1];
    quantities = q.quantities;
    grad = q.grad;
}

bool operator<(const quantity_& q, const double& d) {
    return (q.value[0] < d);
}

bool operator<=(const quantity_& q, const double& d) {
    return (q.value[0] <= d);
}

bool operator<=(const quantity_& q1, const quantity_& q2) {
    return (q1.value[0] <= q2.value[0]);
}

bool operator>(const quantity_& q, const double& d) {
    return (q.value[0] > d);
}

bool operator>=(const quantity_& q, const double& d) {
    return (q.value[0] >= d);
}

void quantity_::check() {
#ifdef DEBUG
    if (isinf(value[0])) {
        warning({"overflow detected in quantity ", name});
    } else if (isnan(value[0])) {
        warning({"Not a number value (NaN) detected in quantity ", name});
    }
#endif
}

//void quantity_::add(const quantity_& q) {
//    quantities = q.quantities;
//
//    for (unsigned int i = 0; i < q.quantities.size(); ++i) {
//        quantity *qi = quantities[i];
//        if (qi == &q)
//            return;
//    }
//    quantities.push_back(q);
//}
//
//void quantity_::add(const quantity_& q1, const quantity_& q2) {
//    quantities = q1.quantities;
//    quantities.insert(q2.quantities.begin(), q2.quantities.end());
//
//    if (quantities.find(q1.name) != quantities.end() &&
//            (q1.type == quantity_::ACROSS || q1.type == quantity_::THROUGH ||
//            q1.type == quantity_::DOUBLE)) {
//        quantities[q1.name] = const_cast<quantity_ *> (&q1);
//    }
//    if (quantities.find(q2.name) != quantities.end() &&
//            (q2.type == quantity_::ACROSS || q2.type == quantity_::THROUGH ||
//            q2.type == quantity_::DOUBLE)) {
//        quantities[q2.name] = const_cast<quantity_ *> (&q2);
//    }
//}

//double quantity_::getGrad(const string& s) const {
//    double retval;
//
//    map<const string, double>::const_iterator ig = grad.find(s);
//
//    if (ig == grad.end()) {
//        //cout << "Warning: " << name << ".grad[" << s << "] does not exist" 
//        //     << endl;
//        retval = 0.0;
//    } else
//        retval = ig->second;
//    return (retval);
//}

void quantity_::init() {
}