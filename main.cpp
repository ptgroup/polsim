#include "simulation.hpp"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    auto sim = Simulation{140.2};

    for (int i = 0; i < 10; i++) {
        cout << sim.take_data() << endl;
        sim.run_for(1);
    }

    return 0;
}
