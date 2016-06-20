#include "util.h"

namespace util {

void report(const string& msg, const int& severity) {
    switch(severity) {
        case WARNING:
            cerr << "WARNING: ";
            break;
        case ERROR:
            cerr << "ERROR: ";
            break;
        case INFO:
        default:
            cerr << "INFO: ";
            break;
    }
    cerr << msg << endl;
}

void report(initializer_list<string> err, const int& severity) {
    switch(severity) {
        case WARNING:
            cerr << "WARNING: ";
            break;
        case ERROR:
            cerr << "ERROR: ";
            break;
        case INFO:
        default:
            cerr << "INFO: ";
            break;
    }
    for (auto& s: err)
        cerr << s;
    cerr << endl;
}

void info(const string& msg) {
    report(msg, INFO);
}

void info(initializer_list<string> err) {
    report(err, INFO);
}

void warning(const string& msg) {
    report(msg, WARNING);
}

void warning(initializer_list<string> err) {
    report(err, WARNING);
}

void error(const string& msg) {
    report(msg, ERROR);
}

void error(initializer_list<string> err) {
    report(err, ERROR);
}

void fatal(const string& msg) {
    report(msg, ERROR);
}

void fatal(initializer_list<string> err) {
    report(err, ERROR);
}

};
