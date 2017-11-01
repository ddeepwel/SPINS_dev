#ifndef BASECASE_IMPL
#define BASECASE_IMPL 1

// append a diagnostic to string for printing into diagnostic file
template <class T> void BaseCase::add_diagnostic(const string str, const T val,
        string & header, string & line) {
    // append to the header
    header.append(str + ", ");
    // append to the line of values
    ostringstream oss;
    oss.precision(17);
    oss << scientific << val;
    line.append(oss.str() + ", ");
}

#endif
