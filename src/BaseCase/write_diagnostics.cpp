#include "../BaseCase.hpp"

// write the diagnostic file
void BaseCase::write_diagnostics(string header, string line,
        int iter, bool restarting) {
    // remove last two elements (comma and space)
    string clean_header = header.substr(0, header.size()-2);
    string clean_line   =   line.substr(0,   line.size()-2);
    // open file
    FILE * diagnos_file = fopen("diagnostics.txt","a");
    assert(diagnos_file);
    // print header
    if (iter == 0 and !restarting) {
        fprintf(diagnos_file,"%s\n",clean_header.c_str());
    }
    // print the line of values
    fprintf(diagnos_file, "%s\n", clean_line.c_str());
    // Close file
    fclose(diagnos_file);
}

// Explicitly instantiate common add_diagnostic templates, for integer and double
// types
template void BaseCase::add_diagnostic<int>(const string str, const int val,
        string & header, string & line);
template void BaseCase::add_diagnostic<double>(const string str, const double val,
        string & header, string & line);
