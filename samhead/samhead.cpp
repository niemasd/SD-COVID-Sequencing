/**
 * Print the first n entries in a given SAM passed via stdin
 */

// include statements
#include <iostream>
#include <sstream>
using namespace std;

// useful constants
#define USAGE_MESSAGE "USAGE: samhead <n> <all/successful/failed>"
#define ERROR_NUM_ARGS "ERROR: Incorrect number of arguments"
#define ERROR_N "ERROR: Invalid value of n"
#define ERROR_SELECTION "ERROR: Invalid selection"

// main function
int main(int argc, char* argv[]) {
    // check number of arguments
    if(argc != 3) {
        cerr << ERROR_NUM_ARGS << endl << USAGE_MESSAGE << endl;
        return 1;
    }

    // parse N (number of SAM entries to output)
    const unsigned long N = strtoul(argv[1], nullptr, 0);
    if(N == 0) {
        cerr << ERROR_N << ": " << argv[1] << endl << USAGE_MESSAGE << endl;
        return 1;
    }

    // parse SAM entry type selection
    bool out_success, out_fail;
    switch(argv[2][0]) {
        case 'a':
        case 'A': out_success = true; out_fail = true; break;
        case 's':
        case 'S': out_success = true; out_fail = false; break;
        case 'f':
        case 'F': out_success = false; out_fail = true; break;
        default: cerr << ERROR_SELECTION << ": " << argv[2] << endl; return 1;
    }

    // parse SAM from stdin
    unsigned long num_success = 0; // number of successful reads seen by samhead
    unsigned long num_fail = 0;    // number of failed reads seen by samhead
    unsigned long num_out = 0;     // number of reads output by samhead
    string line;                   // current line
    string tmp_str_ID;             // temporary holder string for current ID
    string tmp_str_FLAG;           // temporary holder string for bitwise FLAG
    while(num_out < N && getline(cin, line)) {
        // output all header lines
        if(line[0] == '@') {
            cout << line << endl; continue;
        }

        // parse current row
        istringstream ss(line);
        getline(ss, tmp_str_ID, '\t');   // first column (ID)
        getline(ss, tmp_str_FLAG, '\t'); // second column (bitwise FLAG)

        // unmapped read
        if(stoi(tmp_str_FLAG) & 4) {
            ++num_fail;
            if(out_fail) {
                cout << line << endl; ++num_out;
            }
        }

        // mapped read
        else {
            ++num_success;
            if(out_success) {
                cout << line << endl; ++num_out;
            }
        }
    }

    // finish up
    cerr << "Mapped\t" << num_success << endl;
    cerr << "Unmapped\t" << num_fail << endl;

    return 0;
}
