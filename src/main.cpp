#include "polsim/controller.hpp"
#include "polsim/pdp.hpp"
#include "polsim/simulation.hpp"

#include <iostream>

int main()
{
    using namespace polsim;
    using std::cout;
    using std::endl;

    constexpr auto start = 140.1;
    auto controller = StandardController<3>{Pdp{Simulation{start}}};
    auto perfect = PerfectController{Pdp{Simulation{start}}};

    constexpr auto steps = 50;

    cout << "Perfect:" << endl;
    for (auto i = 0; i < 3 * steps; i++)
        cout << perfect.step() << endl;
    cout << "Beam on" << endl;
    perfect.system_ref().beam_on();
    for (auto i = 0; i < 3 * steps; i++)
        cout << perfect.step() << endl;

    cout << "Standard:" << endl;
    for (auto i = 0; i < steps; i++)
        cout << controller.step() << endl;
    cout << "Beam on" << endl;
    controller.system_ref().beam_on();
    for (auto i = 0; i < steps; i++)
        cout << controller.step() << endl;

    return 0;
}
