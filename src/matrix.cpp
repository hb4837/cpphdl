#include "global.h"
#include "node.h"
#include "matrix.h"
#include "util.h"

using namespace util;

matrix::matrix(const int n) {
    m = create(n);
    myerrno = 0;
//    row = 0;
//    col = 0;
    size = n;
    e = 0;
    rhsVector = (vector<double>) 0;
    rhs = new double[n];
    x = new double[n];
    factored = false;
    solved = false;
}

MATRIX *matrix::create(const int n=0) {
    MATRIX *s = spCreate(n, 0, &myerrno);
    // matrixError(err);
    return(s);
}

void matrix::clear() {
    spClear(m);
    clearRhs();
    factored = false;
    solved = false;
    for(int i = 0; i < size; ++i) {
        rhs[i] = 0.0;
        x[i] = 0.0;
    }
}

void matrix::info() {
    int i;

    cout << endl << endl << "RHS: [";

    cout << rhs[0];
    for(i = 1; i < size; ++i)
        cout << ", " << rhs[i];

    cout << "]" << endl << endl;

    /*
     * arguments of spPrint:
     *
     * arg1 = 0: print matrix as input
     * arg1 = 1: print reordered matrix 
     * arg2 = 0: print matrix entries as x's
     * arg2 = 1: print matrix entries as values 
     * arg3 = 0: don't print extra information
     * arg3 = 1: print extra information
     */
    if(!factored)
        spPrint(m, 0, 1, 1, stdout);
    else
        //spPrint(m, 1, 1, 1);
        spPrint(m, 0, 0, 1, stdout);  // without values

    if(solved) {
        cout << endl << endl << "SOLUTION: [";
        cout << x[0];
        for(i = 1; i < size; ++i)
            cout << ", " << x[i];
        cout << "]" << endl << endl << endl;
    }
}

double matrix::findEntry(const int& row, const int& col) {
    e = spFindElement(m, row, col);
    if(e == 0)
        return(0.0);
    else
        return(*e);
}

void matrix::clearEntry(const int& row, const int& col) {
    e = spGetElement(m, row, col);
    *e = 0.0;
}

void matrix::addEntry(const int& row, const int& col) {
    if(row > size)
        size = row;
    if(col > size)
        size = col;
    e = spGetElement(m, row, col);
}

void matrix::addEntry(const int& row, const int& col, const double& val) {
    if (fabs(val) > MACHINE_RESOLUTION) {
        if(row > size)
            size = row;
        if(col > size)
            size = col;
        // cout << "addEntry(" << row << ", " << col << ", " << val << ")" << endl;
        e = spGetElement(m, row+1, col+1);
        spADD_REAL_ELEMENT(e, val);
    }
}

void matrix::setEntry(const int& row, const int& col, const double& val) {
    if (fabs(val) > MACHINE_RESOLUTION) {
        if(row > size)
            size = row;
        if(col > size)
            size = col;
        // cout << "setEntry(" << row << ", " << col << ", " << val << ")" << endl;
        e = spGetElement(m, row, col);
        *e = val;
    }
}

/*
void matrix::addEntry(const node& n1, const node& n2, const double& val) {
    if(n1.ground || n2.ground)
        return;

    int row = n1.index;
    int col = n2.index;

    if(row > size)
        size = row;
    if(col > size)
        size = col;

    e = spGetElement(m, row, col);
    spADD_REAL_ELEMENT(e, val);
}
*/

void matrix::deleteRowAndCol(const int& i) {
    spDeleteRowAndCol(m, i, i);
}

void matrix::addRhsEntry(const int& row) {
    rhs[row] = 0.0;
}

void matrix::addRhsEntry(const int& row, const double& val) {
    if (fabs(val) > MACHINE_RESOLUTION) {
        rhs[row] += val;
        // cout << "addRhsEntry(" << row << ", " << val << ")" << endl;
    }
}

void matrix::setRhsEntry(const int& row, const double& val) {
    if (fabs(val) > MACHINE_RESOLUTION) {
        rhs[row] = val;
        // cout << "setRhsEntry(" << row << ", " << val << ")" << endl;
    }
}

//void matrix::deleteRhsEntry(const int& row) {
//    vector<double>::iterator iv;
//    int i = 1;
//
//    for(iv = rhsVector.begin(); iv != rhsVector.end(); ++iv) {
//        if(i++ == row) {
//            rhsVector.erase(iv);
//            return;
//        }
//    }
//}

void matrix::clearRhs() {
    for (int i = 0; i < size; i++)
        rhs[i] = 0.0;
}

void matrix::factor() { 
    myerrno = spFactor(m); 
    // myerrno = spError(m);
    //if((err = spError(m)) != spOKAY) {
        //cerr << "error in matrix::factor(): " << errString(err) << endl;
        //exit(1);
    //}
    if(myerrno != spOKAY)
        spOrderAndFactor(m, 0, 0, 0, 1); 
    //matrixError(myerrno);
    factored = true;
}

void matrix::solve() {
    int i;

    if(!factored) {
        error({"matrix::solve(): matrix has not been factored yet."});
    }

    for(i = 0; i < size; ++i) {
        x[i] = 0.0;
    }

    spSolve(m, rhs, x);
//  if((myerrno = spError(m)) != spOKAY) {
//      error({"matrix::solve: ", errString(myerrno)});
//  }

    // Check for inf, nan etc.
    for(i = 0; i < size; ++i) {
        if (std::isinf(x[i])) {
            fatal({"x[", to_string(i), "] is infinite."});
        }
        else if (std::isnan(x[i])) {
            fatal({"x[", to_string(i), "] is NaN."});
        }
    }

    solved = true;
}

string matrix::errString(const int& errNumber) {
    string msg("12345678901234567890");

    switch(errNumber) {
        case spOKAY:        msg = "spOKAY"; break;
        case spSMALL_PIVOT: msg = "spSMALL_PIVOT"; break;
        case spZERO_DIAG:   msg = "spZERODIAG"; break;
        case spSINGULAR:    msg = "spSINGULAR"; break;
        case spNO_MEMORY:   msg = "spNO_MEMORY"; break;
        case spPANIC:       msg = "spPANIC"; break;
        default:            msg = "unknown"; break;
    }
    return(msg);
}

