#ifndef _QUANTITY_H
#define _QUANTITY_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "global.h"
#include "model.h"
#include "port.h"
#include "equation.h"

using namespace std;

class across {
};

class through {
};

class quantity_ {
    friend class port;
    friend class model;
    friend class equation;
private:
    int index;
    // Back pointers to the objects containing this quantity
    model *modelPtr;
    equation *equationPtr;
    double valueDelta;
    bool limitFlag;
    bool printFlag;
    bool dynamicFlag;
    bool limited;
protected:
    string name;
    int type;
    double value[2];
    double icValue;
    double dotValue;
    double trValue[4];
    double trDotValue[2];
    vector<quantity_ *> quantities;
    vector<double> grad;
    bool icFlag;
    bool implicitFlag;
    bool dotFlag;
    // Possible types
    static const int UNKNOWN = 0;
    static const int ACROSS = 1;
    static const int THROUGH = 2;
    static const int DOUBLE = 3;
    static const int INTERMEDIATE = 4;
    static const int LOCAL = 5;
public:
    quantity_();
    quantity_(const string&);
    virtual ~quantity_();
    // Getters
    virtual const string& getName() const { return name; }
    inline int getIndex() const { return index; }
    virtual port *getFromPort() const { return 0; }
    virtual port *getToPort() const { return 0; }
    virtual node *getFromNode() const { return 0; }
    virtual node *getToNode() const { return 0; }
    inline bool isAcrossQuantity() const { return (type == quantity_::ACROSS); }
    inline bool isThroughQuantity() const { return (type == quantity_::THROUGH); }
    inline bool isDoubleQuantity() const { return (type == quantity_::DOUBLE); }
    inline double current() const { return value[0]; }
    inline double previous() const { return value[1]; }
    inline double delta() const { return valueDelta; }
    bool isImplicit() const { return (implicitFlag); }
    bool isExplicit() const { return (!implicitFlag); }
    inline double getGrad(const int i) const { return grad[i]; }
    inline double getValue(const int i=0) const { return value[i]; }
    inline equation *getEquation() const { return equationPtr; }
    inline model *getModel() const { return modelPtr; }
    // Setters
    void setIndex(int i) { index = i; }
    virtual void setFromPort(port *p) {};
    virtual void setToPort(port *p) {};
    friend ostream& operator<<(ostream&, const quantity_&);
    inline void setDynamicFlag() { dynamicFlag = true; }
//    void add(const quantity_&);
//    void add(const quantity_&, const quantity_&);
    void setExplicit() { implicitFlag = false; }
    void setImplicit() { implicitFlag = true; }
    void setValue(double v) { value[0] = v; }
    void pushValue(double v) { value[1] = value[0]; value[0] = v; }
    void pushAddValue(double v) { value[1] = value[0]; value[0] += v; }
    void init();
    void reset();
    void check();
    // Operators
//    friend quantity_ pow(const quantity_&, const double&);
//    friend quantity_ exp(const quantity_&);
//    friend quantity_ limexp(const quantity_&);
//    friend quantity_ log(const quantity_&);
//    friend quantity_ sqrt(const quantity_&);
//    friend quantity_ sin(const quantity_&);
    friend quantity_ operator-(const quantity_&);
    friend quantity_ operator+(const quantity_&, const quantity_&);
    friend quantity_ operator+(const quantity_&, const double&);
    friend quantity_ operator+(const double&, const quantity_&);
    friend quantity_ operator-(const quantity_&, const quantity_&);
    friend quantity_ operator-(const quantity_&, const double&);
    friend quantity_ operator-(const double&, const quantity_&);
    friend quantity_ operator*(const quantity_&, const quantity_&);
    friend quantity_ operator*(const quantity_&, const double&);
    friend quantity_ operator*(const double&, const quantity_&);
    friend quantity_ operator/(const quantity_&, const quantity_&);
    friend quantity_ operator/(const quantity_&, const double&);
    friend quantity_ operator/(const double&, const quantity_&);
    friend void operator==(const quantity_&, const quantity_&);
    friend void operator==(const quantity_&, const double&);
    void operator=(const quantity_&);
    friend bool operator<(const quantity_&, const double&);
    friend bool operator<=(const quantity_&, const double&);
    friend bool operator<=(const quantity_&, const quantity_&);
    friend bool operator>(const quantity_&, const double&);
    friend bool operator>=(const quantity_&, const double&);
    quantity_ dot();
    friend double double_cast(const quantity_&);
};

template<class T> class quantity : public quantity_ {
private:
    port *from;
    port *to;
public:
    quantity();
    quantity(const string&);
    quantity(const string&, port&);
    quantity(const string&, port&, port&);
    quantity(quantity_);
    // Getters
    inline node *getFromNode() const { return from->nodePtr; }
    inline node *getToNode() const { return to->nodePtr; }
    inline port *getFromPort() const { return from; }
    inline port *getToPort() const { return to; }
    // Setters
    void setFromPort(port *p) { from = p; }
    void setToPort(port *p) { to = p; }
    void startvalue(const double&);
    // Operators
    void operator=(const double&);
    //void operator =(const quantity_&);
};

#endif // _QUANTITY_H
