#ifndef _UTIL_H
#define _UTIL_H

#include <iostream>
#include <string>

using namespace std;

namespace util {

static const int INFO = 0;
static const int WARNING = 1;
static const int ERROR = 2;

void report(const string&, const int& severity);
void report(initializer_list<string>, const int& severity);

void error(const string&);
void error(initializer_list<string>);

void warning(const string&);
void warning(initializer_list<string>);

void fatal(const string&);
void fatal(initializer_list<string>);

} /* namespace util */

#endif // _UTIL_H
