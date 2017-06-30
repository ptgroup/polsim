#include "controller.hpp"
#include "pdp.hpp"
#include "simulation.hpp"

#include <iostream>

int main(int argc, char **argv)
{
    using std::cout;
    using std::endl;

    auto controller = StandardController<4>{Pdp{Simulation{140.0}}};
    auto perfect = PerfectController{Pdp{Simulation{140.0}}};

    cout << "Perfect:" << endl;
    for (auto i = 0; i < 40; i++)
        cout << perfect.step() << endl;

    cout << "Standard:" << endl;
    for (auto i = 0; i < 10; i++)
        cout << controller.step() << endl;

    return 0;
}
