#ifndef PORT_H
#define PORT_H

#include <iostream>
#include <iomanip>
#include <strstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include "node.h"

using namespace std;

class model;

class port {
    friend class node;
    friend class model;
    template<class T> friend class quantity;
private:
    string name;
    node *nodePtr;
public:
    port();
    port(string);
    port(node&);
    virtual ~port() { }
    void add(port *);
    inline node *getNode() const { return nodePtr; }
    void operator() (node&);
    friend ostream& operator <<(ostream&, const port&);
};

#endif /* PORT_H */

