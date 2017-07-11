/**
 * @file StandardController.tpp
 * @brief Template methods for StandardController2.
 *
 * @author Ian Johnson
 * @date July 11, 2017
 */
#include "controller.hpp"

template <unsigned n_points>
StandardController2<n_points>::StandardController2(Pdp pdp, double step_size,
                                                   bool seek_positive)
    : Controller(pdp, seek_positive), step_size(step_size)
{
}

template <unsigned n_points>
Data StandardController2<n_points>::step()
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

	// Compare the rate to the last one and adjust the step size
	// accordingly.
	if (avg_rate < this->last_rates.back()) {
		this->step_size *= STEP_SIZE_REDUCE;
		this->direction *= -1.0;
		// We shouldn't go below a certain step size.
		if (this->step_size < MIN_STEP_SIZE)
			this->step_size = MIN_STEP_SIZE;
	}
	// We may need to use the correctional factor if the collected rate is
	// less than the minimum of previous rates times some fraction.  The
	// BAD_FRACTION means that it doesn't necessarily have to be below the
	// minimum rate, but only below a certain fraction of the minimum rate
	// which is considered "bad enough"; alternatively, BAD_FRACTION could
	// be made greater than 1, in which case only rates which are even
	// worse than the minimum would be considered bad.
	if (this->last_rates.size() == N_RATES &&
	    avg_rate < this->minimum_rate() * BAD_FRACTION) {
		std::cout << "correcting" << std::endl;
		this->step_size *= STEP_SIZE_CORRECT;
		// Don't correct too much; clear the list of rates to prevent
		// overcorrection.
		this->last_rates.clear();
	}

	// Move the "motor" and store the last rate.
	this->pdp.set_freq(last_data.freq + this->direction * this->step_size);
	this->add_rate(avg_rate);

	return last_data;
}

template <unsigned n_points>
void StandardController2<n_points>::add_rate(double rate)
{
	// Add the new rate, getting rid of the oldest one if necessary.
	this->last_rates.push_back(rate);
	if (this->last_rates.size() > N_RATES)
		this->last_rates.pop_front();
}

template <unsigned n_points>
double StandardController2<n_points>::minimum_rate()
{
	auto min_rate = 0.0;
	for (auto rate : this->last_rates)
		if (rate < min_rate)
			min_rate = rate;
	return min_rate;
}
