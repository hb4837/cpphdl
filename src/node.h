#ifndef _NODE_H
#define _NODE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <strstream>

using namespace std;

#define GND(x) node gnd()
#define NODE(x) node x(#x)

class node {
    friend class quantity_;
    friend class equation;
    friend class model;
    friend class matrix;
    friend class port;
private:
    string name;
    int index;
    double value[2];
    bool icFlag;
    double icValue;
    double trValue[4];
    double dotValue;
    double trDotValue[2];
    double valueDelta;
    int ncon;
    bool ground;
    bool printFlag;
public:
    node();
    node(string);
    node(string, double);
    // Getters
    inline int getIndex() const { return index; }
    inline bool isGround() const { return ground; }
    inline double getValue(const int& i=0) const { return value[i]; }
    // Setters
    inline void setIndex(const int i) { index = i; }
    inline void setValue(const double d) { value[0] = d; }
    inline void pushValue(const double d) { value[1] = value[0]; value[0] = d; }
    inline void pushAddValue(const double d) { value[1] = value[0]; value[0] += d; }
    void reset();
    friend ostream& operator <<(ostream&, const node&);
    const string& getName() const { return(name); }
    void init();
    void addConnection();
};
#endif // _NODE_H
