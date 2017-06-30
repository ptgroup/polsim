#include "controller.hpp"
#include "pdp.hpp"
#include "simulation.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    auto controller = StandardController<4>{Pdp{Simulation{140.2}}};

    for (int i = 0; i < 10; i++) {
        std::cout << controller.step() << std::endl;
    }

    return 0;
}
