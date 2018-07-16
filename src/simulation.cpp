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

void System::beam_on(double current) { this->beam_current = current; }

void System::beam_off() { this->beam_current = 0; }

std::ostream &operator<<(std::ostream &out, const Data &data)
{
    // Output in the standard CSV format
    return out << (int64_t)data.t << "," << data.pn << "," << data.freq;
}

Simulation::Simulation(double freq)
{
    this->set_freq(freq);
    // Make sure pe0 gets initialized
    this->set_temperature(this->system.temperature);
    this->update_transition_rates();
    this->rng.seed(std::random_device{}());
}

Data Simulation::take_data() const
{
    return {this->t, this->pn,          this->pe,  this->freq,
            this->c, this->temperature, this->dose};
}

System &Simulation::system_ref() { return this->system; }

void Simulation::set_freq(double frequency)
{
    this->freq = frequency;
    // Make sure transition rates get updated
    this->update_transition_rates();
}

void Simulation::run_for(double t, double step)
{
    const double end_t = this->t + t;
    while (this->t < end_t)
        this->time_step(step);
    if (t > 100.0)
        std::cout << fit_params.a << std::endl;
}

void Simulation::anneal(double t, double temperature)
{
    // Reset phi (remove negative effects of irradiation)
    this->phi = 0;
    // Maybe change T1n?
    this->t1n *= 0.8;

    auto temp_tmp = this->system.temperature;
    this->system.set_temperature(temperature);
    this->run_for(t);
    this->system.set_temperature(temp_tmp);
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
    this->pe0 = -tanh(2 / temperature);
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
        this->system.temperature + 0.25 * this->system.beam_current / 100;

    // Increase phi linearly.
    // const double k_phi = (this->system.beam_current - 30) / 1e8;
    const double k_phi = 0;

    // Calculate convenience constants
    const double A = -this->t1e / this->t1n -
                     (this->c / 2) * (this->alpha + this->beta) - this->phi;
    const double B = (this->c / 2) * (this->alpha - this->beta);
    const double C = (this->alpha - this->beta) / 2;
    const double D = -1 - (this->alpha + this->beta) / 2;

    // Calculate rates
    const double pn_prime = (A * this->pn_raw + B * this->pe) / this->t1e;
    const double pe_prime =
        (C * this->pn_raw + D * this->pe + this->pe0) / this->t1e;

    // Update pn and pe using Euler's method
    this->pn_raw += pn_prime * t;
    this->pe += pe_prime * t;
    // Update temperature and phi
    this->set_temperature(this->temperature +
                          t * k_temp * (temp_ss - this->temperature));
    this->phi += t * k_phi;
    if (this->phi < 0)
        this->phi = 0;
    // Update C
    this->c += IRRADIATION_FACTOR * this->system.beam_current * t;

    // Update dose, along with the fit parameters M1 and M2, which depend on
    // it.
    // Recall the exponential change of M1 and M2 described in the
    // documentation.
    const double delta_dose =
        BEAM_FRACTION * (this->system.beam_current * 1e-9 / ELEM_CHARGE) * t;
    this->fit_params.m1 +=
        FIT_M1_COEFF * FIT_M1_RATE * delta_dose * exp(FIT_M1_RATE * this->dose);
    this->fit_params.m2 +=
        FIT_M2_COEFF * FIT_M2_RATE * delta_dose * exp(FIT_M2_RATE * this->dose);
    this->fit_params.a += FIT_A_RATE * delta_dose * this->fit_params.a;
    this->dose += (this->system.beam_current * 1e-9 / ELEM_CHARGE) * t;

    // Calculate new transition rates
    this->update_transition_rates();
    // Update time
    this->t += t;
    // Update "noisy pn"
    this->pn = this->pn_noisy();
}

void Simulation::update_transition_rates()
{
    auto rates = this->calc_transition_rates(this->freq);
    this->alpha = rates.first;
    this->beta = rates.second;
}

std::pair<double, double> Simulation::calc_transition_rates(double freq) const
{
    // Calculate from the Gaussian distributions
    const double scale =
        this->fit_params.a / (sqrt(2 * PI) * this->fit_params.s);
    const double diff1 = freq - this->fit_params.m1;
    const double diff2 = freq - this->fit_params.m2;
    const double exp1 =
        exp(-diff1 * diff1 / (2 * this->fit_params.s * this->fit_params.s));
    const double exp2 =
        exp(-diff2 * diff2 / (2 * this->fit_params.s * this->fit_params.s));

    return std::make_pair(scale * exp2, scale * exp1);
}

double Simulation::pn_noisy()
{
    auto dist = std::uniform_real_distribution<double>{-1, 1};
    const double THERMAL_NOISE = THERMAL_RANDOMNESS * dist(this->rng);
    const double UNIFORM_NOISE = BASE_RANDOMNESS * dist(this->rng);

    return this->pn_raw * (1 + THERMAL_NOISE) + UNIFORM_NOISE;
}

double Simulation::steady_state(double freq) const
{
    auto rates = this->calc_transition_rates(freq);
    auto alpha = rates.first;
    auto beta = rates.second;
    // See Jeffries-optimal-frequency.nb for details
    return this->c * this->t1n * (beta - alpha) /
           (this->t1e * (2 + alpha + beta) +
            this->c * this->t1n * (alpha + beta + 2 * alpha * beta));
}
} // namespace polsim
