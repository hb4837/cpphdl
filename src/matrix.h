
#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
//#include "../sparse/spDefs.h"
#define spINSIDE_SPARSE
#include "../sparse/spConfig.h"
#undef spINSIDE_SPARSE
#include "../sparse/spMatrix.h"

typedef char MATRIX;

class matrix {
private:
    int myerrno;
    MATRIX *m;
//    int row, col;
    int size;
    double *e;
    vector<double> rhsVector;
    double *rhs;
    double *x;
    bool factored;
    bool solved;
public:
    matrix(const int);
    MATRIX *create(const int);
    void clear(); 
    void info();
    double findEntry(const int&, const int&);
    void clearEntry(const int&, const int&);
    void addEntry(const int&, const int&);
    void addEntry(const int&, const int&, const double&);
    void addEntry(const node&, const node&, const double&);
    void setEntry(const int&, const int&, const double&);
    void deleteRowAndCol(const int&);
    void addRhsEntry(const int&);
    void setRhsEntry(const int&, const double&);
    void addRhsEntry(const int&, const double&);
    void deleteRhsEntry(const int&);
    void factor();
    void solve();
    inline double getSolution(const int&i) const { return(x[i]); }
    void clearRhs();
    string errString(const int&);
    int getSize() const { return(size); }
};
#endif // _MATRIX_H




