#include <string>
#include <stdio.h>
#include <vector>
using namespace std;

std::vector<std::string> split_string(const string& str, const char& ch) 
{
    string next("");
    vector<string> result;

    for (string::const_iterator it = str.begin(); it != str.end(); it++) {
        // if the delimiter character is hit
        if (*it == ch) {
            // if some characters are accumulated
            if (!next.empty()) {
                // they are added to the result vector
                result.push_back(next);
                next.clear();
            }
        } else {
            // the next character is accumulated into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(next);
    return result;
}
