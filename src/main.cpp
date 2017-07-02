#include "polsim/controller.hpp"
#include "polsim/pdp.hpp"
#include "polsim/simulation.hpp"

#include <iostream>

int main()
{
	using namespace std;
	using namespace polsim;

	auto controller = StandardController<4>{Pdp{Simulation{140.0}}};
	controller.system_ref().beam_on();
	auto perfect = PerfectController{Pdp{Simulation{140.0}}};
	controller.system_ref().beam_on();

	cout << "Perfect:" << endl;
	for (auto i = 0; i < 40; i++)
		cout << perfect.step() << endl;

	cout << "Standard:" << endl;
	for (auto i = 0; i < 10; i++)
		cout << controller.step() << endl;

	return 0;
}
