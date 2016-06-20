#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <iostream>
#include <iomanip>
#include <strstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

using namespace std;

class parameter_ {
    friend class kernel;
protected:
    string name;
public:
    parameter_();
    parameter_(const string);

    virtual ~parameter_() {
    }
    void insert(parameter_ *);
    friend ostream& operator<<(ostream&, const parameter_&);
};

template<class T> class parameter : public parameter_ {
private:
    bool defFlag;
    T value;
public:
    parameter() : parameter_() {
        defFlag = false;
    }

    parameter(const string s) : parameter_(s) {
        defFlag = false;
    }

    parameter(string s, const T& v) : parameter_(s) {
        value = v;
        defFlag = true;
    }
    
    inline T getValue() const { return value; }
    inline bool defined() const { return defFlag; }
    inline void operator()(const T& v) { value = v; defFlag = true; }
    void operator=(const T& v) { value = v; }
    inline bool operator==(const T& v) { return (value == v); }
    inline bool operator!=(const T& s) { return (!(value == s)); }
    inline bool operator<=(const T& s) { return ((value <= s)); }
    inline bool operator>=(const T& s) { return (!(value >= s)); }
    inline operator double() { return value; }

    ostream& operator<<(ostream& os) {
        os << name << " " << value;
        return os;
    }
};

#endif // _PARAMETER_H
