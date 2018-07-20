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
    : NPointController<n_points>(pdp, step_size, seek_positive)
{
}

template <unsigned n_points>
typename NPointController<n_points>::Decision
StandardController<n_points>::make_decision(
    const std::array<Data, n_points> &data)
{
    return seek_pol ? make_decision_polarization(data)
                    : make_decision_rate(data);
}

template <unsigned n_points>
typename NPointController<n_points>::Decision
StandardController<n_points>::make_decision_rate(
    const std::array<Data, n_points> &data)
{
    double k_avg = 0.0;

    for (unsigned i = 2; i < data.size(); i++) {
        const auto t0 = data[i - 2].t, t1 = data[i - 1].t, t2 = data[i].t;
        const auto p0 = data[i - 2].pn, p1 = data[i - 1].pn, p2 = data[i].pn;
        const auto c0 = p0 / ((t0 - t1) * (t0 - t2)),
                   c1 = p1 / ((t1 - t0) * (t1 - t2)),
                   c2 = p2 / ((t2 - t0) * (t2 - t1));
        const auto pn_prime = [=](double t) {
            return c0 * (2 * t - t1 - t2) + c1 * (2 * t - t0 - t2) +
                   c2 * (2 * t - t0 - t1);
        };

        k_avg += pn_prime(t1);
    }

    k_avg /= data.size() - 2;

    // We negate the rates if we're seeking negative to make the calculations
    // easier (no need for separate cases).
    const auto k_before = this->seek_positive ? last_k : -last_k;
    const auto k_after = this->seek_positive ? k_avg : -k_avg;
    using Decision = typename NPointController<n_points>::Decision;
    Decision decision;
    if (k_before <= 0.0 || k_after <= 0.0) {
        decision = k_after >= k_before ? Decision::KEEP_DIRECTION
                                       : Decision::SWITCH_DIRECTION;
    } else {
        decision = fabs(k_after / k_before) >= GOOD_RATIO
                       ? Decision::KEEP_DIRECTION
                       : Decision::SWITCH_DIRECTION;
    }

    last_k = k_avg;
    return decision;
}

template <unsigned n_points>
typename NPointController<n_points>::Decision
StandardController<n_points>::make_decision_polarization(
    const std::array<Data, n_points> &data)
{
    const auto avg = std::accumulate(data.begin(), data.end(), 0,
                                     [&](auto a, auto b) { return a + b.pn; }) /
                     data.size();
    using Decision = typename NPointController<n_points>::Decision;
    const auto decision =
        avg >= last_pol ? Decision::KEEP_DIRECTION : Decision::SWITCH_DIRECTION;
    last_pol = avg;

    return decision;
}

template <unsigned n_points>
Data StandardController<n_points>::step()
{
    const auto data = this->NPointController<n_points>::step();
    if (data.t > ALGO_SWITCH_TIME) {
        seek_pol = true;
    }
    return data;
}
