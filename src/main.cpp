#include "polsim/controller.hpp"
#include "polsim/pdp.hpp"
#include "polsim/simulation.hpp"

#include <iostream>

int main()
{
	using namespace std;
	using namespace polsim;

	auto controller = StandardController<3>{Pdp{Simulation{140.0}}};
	auto perfect = PerfectController{Pdp{Simulation{140.0}}};
	auto experimental = StandardController2<3>{Pdp{Simulation{140.0}}};
	controller.system_ref().beam_on();
	perfect.system_ref().beam_on();
	experimental.system_ref().beam_on();

	constexpr auto steps = 50;

	cout << "Perfect:" << endl;
	for (auto i = 0; i < 3 * steps; i++)
		cout << perfect.step() << endl;

	cout << "Standard:" << endl;
	for (auto i = 0; i < steps; i++)
		cout << controller.step() << endl;

	cout << "Experimental:" << endl;
	for (auto i = 0; i < steps; i++)
		cout << experimental.step() << endl;

	return 0;
}
