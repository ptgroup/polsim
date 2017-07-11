/**
 * @file StandardController.tpp
 * @brief Template methods for StandardController.
 *
 * @author Ian Johnson
 * @date June 30, 2017
 */
#include "controller.hpp"

template <unsigned n_points>
StandardController<n_points>::StandardController(Pdp pdp, double step_size,
                                                 bool seek_positive)
    : Controller(pdp, seek_positive), step_size(step_size)
{
}

template <unsigned n_points>
Data StandardController<n_points>::step()
{
	// Collect the specified number of data points and average the rates
	auto avg_rate = 0.0;
	auto last_data = this->pdp.take_data();
	for (auto i = 1U; i < n_points; i++) {
		auto new_data = this->pdp.take_data();
		// Calculate the rate
		avg_rate +=
		    (new_data.pn - last_data.pn) / (new_data.t - last_data.t);
		last_data = new_data;
	}
	avg_rate /= n_points;

	// Compare the rate to the last one and move accordingly
	if (avg_rate < last_rate) {
		this->step_size *= STEP_SIZE_REDUCE;
		this->direction *= -1.0;
		// We shouldn't go below a certain step size.
		if (this->step_size < MIN_STEP_SIZE)
			this->step_size = MIN_STEP_SIZE;
	}
	this->pdp.set_freq(last_data.freq + this->direction * this->step_size);
	last_rate = avg_rate;

	return last_data;
}
