#include "polsim/simulation.hpp"

#include <cmath>
#include <iostream>

#include <boost/math/tools/minima.hpp>

namespace polsim
{
void System::set_temperature(double temperature)
{
    this->temperature = temperature;
}

void System::beam_on(double current) { beam_current = current; }

void System::beam_off() { beam_current = 0; }

std::ostream &operator<<(std::ostream &out, const Data &data)
{
    // Output in the standard CSV format
    return out << (int64_t)data.t << "," << data.pn << "," << data.freq;
}

Simulation::Simulation(double freq, FitParameters fit_params)
{
    set_freq(freq);
    // Make sure pe0 gets initialized
    set_temperature(system.temperature);
    set_fit_params(fit_params);
    rng.seed(std::random_device{}());
}

Data Simulation::take_data() const
{
    return {t, pn, pe, freq, fit_params.c, temperature, dose};
}

System &Simulation::system_ref() { return system; }

void Simulation::set_fit_params(FitParameters params)
{
    fit_params = params;
    update_transition_rates();
}

void Simulation::set_freq(double frequency)
{
    freq = frequency;
    // Make sure transition rates get updated
    update_transition_rates();
}

void Simulation::run_for(double t, double step)
{
    const double end_t = this->t + t;
    while (this->t < end_t)
        time_step(step);
    if (t > 100.0)
        std::cout << fit_params.a << " " << std::scientific << dose
                  << std::fixed << std::endl;
}

void Simulation::anneal(double t, double temperature)
{
    // Reset phi (remove negative effects of irradiation)
    phi = 0;
    // Maybe change T1n?
    fit_params.t1n *= 0.8;

    auto temp_tmp = system.temperature;
    system.set_temperature(temperature);
    run_for(t);
    system.set_temperature(temp_tmp);
}

double Simulation::find_optimal_freq(bool negative) const
{
    using boost::math::tools::brent_find_minima;

    // There's only a function to find the minimum, so we'll need to adjust
    // our
    // lambda for negative versus positive polarization
    auto min = negative
                   ? brent_find_minima(
                         [&](double freq) { return this->steady_state(freq); },
                         139.0, 141.0, 20)
                   : brent_find_minima(
                         [&](double freq) { return -this->steady_state(freq); },
                         139.0, 141.0, 20);

    return min.first;
}

void Simulation::set_temperature(double temperature)
{
    // The parameter Pe0 depends on temperature
    pe0 = -tanh(2 / temperature);
    this->temperature = temperature;
}

void Simulation::time_step(double t)
{
    // Parameters for temperature change (exponential growth/decay)
    // temp_ss = steady-state temperature
    // k_temp = rate of exponential increase
    // If we're annealing, we shouldn't allow the temperature to change
    // (assume
    // anneals occur at constant temperature)
    const double k_temp = 1;
    const double temp_ss =
        system.temperature + 0.25 * system.beam_current / 100;

    // Increase phi linearly.
    const double k_phi = 5e-9 * (system.beam_current / 100.0);

    // Calculate convenience constants
    const double A = -fit_params.t1e / fit_params.t1n -
                     (fit_params.c / 2) * (alpha + beta) - phi;
    const double B = (fit_params.c / 2) * (alpha - beta);
    const double C = (alpha - beta) / 2;
    const double D = -1 - (alpha + beta) / 2;

    // Calculate rates
    const double pn_prime = (A * pn_raw + B * pe) / fit_params.t1e;
    const double pe_prime = (C * pn_raw + D * pe + pe0) / fit_params.t1e;

    // Update pn and pe using Euler's method
    pn_raw += pn_prime * t;
    pe += pe_prime * t;
    // Update temperature and phi
    set_temperature(temperature + t * k_temp * (temp_ss - temperature));
    phi += t * k_phi;
    // Update C
    fit_params.c += IRRADIATION_FACTOR * system.beam_current * t;

    // Update dose, along with the fit parameters M1 and M2, which depend on
    // it.
    // Recall the exponential change of M1 and M2 described in the
    // documentation.
    const double delta_dose = (system.beam_current * 1e-9 / ELEM_CHARGE) * t;
    fit_params.m1 +=
        FIT_M1_COEFF * FIT_M1_RATE * delta_dose * exp(FIT_M1_RATE * dose);
    fit_params.m2 +=
        FIT_M2_COEFF * FIT_M2_RATE * delta_dose * exp(FIT_M2_RATE * dose);
    fit_params.a += FIT_A_RATE * delta_dose * fit_params.a;
    dose += delta_dose;

    // Calculate new transition rates
    update_transition_rates();
    // Update time
    this->t += t;
    // Update "noisy pn"
    pn = pn_noisy();
}

void Simulation::update_transition_rates()
{
    auto rates = calc_transition_rates(freq);
    alpha = rates.first;
    beta = rates.second;
}

std::pair<double, double> Simulation::calc_transition_rates(double freq) const
{
    // Calculate from the Gaussian distributions
    const double scale = fit_params.a / (sqrt(2 * PI) * fit_params.s);
    const double diff1 = freq - fit_params.m1;
    const double diff2 = freq - fit_params.m2;
    const double exp1 = exp(-diff1 * diff1 / (2 * fit_params.s * fit_params.s));
    const double exp2 = exp(-diff2 * diff2 / (2 * fit_params.s * fit_params.s));

    return std::make_pair(scale * exp2, scale * exp1);
}

double Simulation::pn_noisy()
{
    auto dist = std::uniform_real_distribution<double>{-1, 1};
    const double THERMAL_NOISE = THERMAL_RANDOMNESS * dist(rng);
    const double UNIFORM_NOISE = BASE_RANDOMNESS * dist(rng);

    return pn_raw * (1 + THERMAL_NOISE) + UNIFORM_NOISE;
}

double Simulation::steady_state(double freq) const
{
    auto rates = calc_transition_rates(freq);
    auto alpha = rates.first;
    auto beta = rates.second;
    // See Jeffries-optimal-frequency.nb for details
    return fit_params.c * fit_params.t1n * (beta - alpha) /
           (fit_params.t1e * (2 + alpha + beta) +
            fit_params.c * fit_params.t1n * (alpha + beta + 2 * alpha * beta));
}
} // namespace polsim
