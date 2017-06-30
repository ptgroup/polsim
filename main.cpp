#include "pdp.hpp"
#include "simulation.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    auto pdp = Pdp{Simulation{140.2}};

    for (int i = 0; i < 10; i++) {
        std::cout << pdp.take_data() << std::endl;
    }

    return 0;
}
