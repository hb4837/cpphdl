#ifndef _EQUATION_H
#define _EQUATION_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>

using namespace std;

class quantity_;
class node;
class model;

class equation {
    friend class model;
private:
    string name;
    int type;
    int index;
    model *modelPtr;
    double f;
    vector<quantity_ *> quantities;
    map<quantity_ *, double> j;
    // Possible types
    static const int NONE = 0;
    static const int NODAL = 1;
    static const int BRANCH = 2;
public:
    equation();
    equation(quantity_&);
    equation(quantity_&, quantity_&);
    ~equation();
    // Setters
    void insert(quantity_ *);
    inline void setIndex(const int i) { index = i; }
    inline void setF(double v) { f = v; }
    inline void setJ(quantity_& q, double v) { j[&q] = v; }
    // Getters
    inline const string& getName() const { return(name); }
    inline const vector<quantity_ *>& getQuantities() const { return quantities; }
    inline double getF() const { return f; }
    inline double getJ(quantity_& q)  { return j[&q]; }
    inline int getIndex() const { return index; }
    // Misc
    void setup();
    void init();
    void load();
    friend ostream& operator <<(ostream&, const equation&);
    // Operators
    void operator =(const quantity_&);
    friend void operator ==(const quantity_&, const quantity_&);
    friend void operator ==(const quantity_&, const double&);
};
#endif // _EQUATION_H

