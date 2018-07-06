/**
 * @file ExperimentalController.tpp
 * @brief Template methods for ExperimentalController.
 *
 * @author Ian Johnson
 * @date June 30, 2017
 */
#include "controller.hpp"

template <unsigned n_points>
ExperimentalController<n_points>::ExperimentalController(Pdp pdp,
                                                         double step_size,
                                                         bool seek_positive)
    : Controller(pdp, seek_positive), step_size(step_size)
{
}

template <unsigned n_points>
Data ExperimentalController<n_points>::step()
{
	// Collect the specified number of data points.
	auto data = std::vector<Data>{};
	for (auto i = 0U; i < n_points; i++) {
		data.push_back(pdp.take_data());
	}
	const auto last = data.back();
	const auto first_t = data.front().t;
	std::for_each(data.begin(), data.end(),
	              [&](auto &point) { point.t -= first_t; });

	const auto t0 = data[0].t, t1 = data[1].t, t2 = data[2].t;
	const auto p0 = data[0].pn, p1 = data[1].pn, p2 = data[2].pn;
	const auto c0 = p0 / ((t0 - t1) * (t0 - t2)),
	           c1 = p1 / ((t1 - t0) * (t1 - t2)),
	           c2 = p2 / ((t2 - t0) * (t2 - t1));
	const auto pn = [=](double t) {
		return c0 * (t - t1) * (t - t2) + c1 * (t - t0) * (t - t2) +
		       c2 * (t - t0) * (t - t1);
	};
	const auto pn_prime = [=](double t) {
		return c0 * (2 * t - t1 - t2) + c1 * (2 * t - t0 - t2) +
		       c2 * (2 * t - t0 - t1);
	};
	const auto pn_2prime = [=](double t) {
		return 2 * c0 + 2 * c1 + 2 * c2;
	};

	const auto k = pn_prime(t1);

	// Compare the rate to the last one and move accordingly
	if (k / last_k < 0.9) {
		this->step_size *= STEP_SIZE_REDUCE;
		this->direction *= -1.0;
		// We shouldn't go below a certain step size.
		if (this->step_size < MIN_STEP_SIZE)
			this->step_size = MIN_STEP_SIZE;
	}
	this->pdp.set_freq(data.back().freq +
	                   this->direction * this->step_size);
	last_k = k;

	return last;
}
