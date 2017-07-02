#include "polsim/controller.hpp"

namespace polsim
{
Controller::Controller(Pdp pdp, bool seek_positive)
    : pdp(pdp), seek_positive(seek_positive)
{
}

System &Controller::system_ref() { return this->pdp.system_ref(); }

PerfectController::PerfectController(Pdp pdp, bool seek_positive)
    : Controller(pdp, seek_positive)
{
}

Data PerfectController::step()
{
	// Take data and then move to the new optimal frequency
	auto data = this->pdp.take_data();
	this->pdp.set_freq(
	    this->pdp.sim.find_optimal_freq(!this->seek_positive));
	return data;
}

ExperimentalController::ExperimentalController(Pdp pdp, bool seek_positive)
    : Controller(pdp, seek_positive)
{
}

Data ExperimentalController::step()
{
	//
}
}
