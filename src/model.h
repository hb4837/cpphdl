#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include "equation.h"
#include "quantity.h"
#include "port.h"
#include "parameter.h"

#define MODEL(x) class x : public model

using namespace std;

class model {
    friend class node;
    friend class quantity_;
    friend class equation;
    friend class parameter_;
    friend class port;
private:
    list<node *> nodes;
    list<parameter_ *> params;
    list<quantity_ *> quantities;
    list<equation *> equationsList;
    list<port *> ports;
    static list<model *> instances;
    static model *currentModel;
    // Constants
    static const int INIT = 0;
    static const int EQUATIONS = 1;
protected:
    string name;
public:
    model(const string);
    model(const char *);
    virtual ~model() { }
    // Getters
    inline const string& getName() const { return name; }
    inline const list<port *>& getPorts() const { return ports; }
    inline const list<node *>& getNodes() const { return nodes; }
    inline const list<quantity_ *>& getQuantities() const { return quantities; }
    inline const list<equation *>& getEquations() const { return equationsList; }
    inline static const list<model *>& getInstances() { return instances; }
    // Setters
    // Misc
    void init();
    virtual void parameters();
    virtual void equations();
    virtual void load();
    virtual void limit();
    virtual void sweep(const double&) { }
    // Operators
    friend ostream& operator <<(ostream&, const model&);
    // Static methods
    static void add(node *);
    static void add(parameter_ *);
    static void add(quantity_ *);
    static void add(equation *);
    static void add(port *);
    static void printInstances();
};
#endif // _MODEL_H
