#include <iostream>
#include "include/bigint.hpp"

using MPI = zeno::MultiPrecisionInteger;

int main() {

    MPI a(123456789);
    MPI b; std::cin >> b; // e.g. b = 98765432143243204378204578

    MPI c("-123456789123456789123456789");

    // Basic Arithmetic Operations
    MPI addition       = a + b;
    MPI subtraction    = a - b;
    MPI multiplication = a * b;
    MPI division       = c / a; 
    MPI remainder      = c % a;

    // Increment and Decrement
    ++a, a--;

    // Comparison Operations and outputs
    std::cout << "b + c: " << b + c << std::endl;
    std::cout << "Is a equal to b? " << (a == b ? "Yes" : "No") << std::endl;
    std::cout << "Is a less than b? " << (a < b ? "Yes" : "No") << std::endl;

    return 0;
}