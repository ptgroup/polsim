#include "simulation.hpp"

#include <cmath>
#include <iostream>

#include <boost/math/tools/minima.hpp>

std::ostream &operator<<(std::ostream &out, const Data &data)
{
    // Output in the standard CSV format
    return out << (int64_t)data.t << "," << data.pn << "," << data.freq;
}

Simulation::Simulation(double freq)
{
    this->set_freq(freq);
    // Make sure pe0 gets initialized
    this->set_temperature(this->system_temperature);
    this->update_transition_rates();
    this->rng.seed(std::random_device{}());
}

Data Simulation::take_data()
{
    return {this->t, this->pn,          this->pe,  this->freq,
            this->c, this->temperature, this->dose};
}

void Simulation::set_freq(double frequency)
{
    this->freq = frequency;
    // Make sure transition rates get updated
    this->update_transition_rates();
}

void Simulation::set_system_temperature(double temperature)
{
    this->system_temperature = temperature;
}

void Simulation::beam_on(double current) { this->beam_current = current; }

void Simulation::beam_off() { this->beam_current = 0; }

void Simulation::run_for(double t, double step)
{
    const double END_T = this->t + t;
    while (this->t < END_T)
        this->time_step(step);
}

void Simulation::anneal(double t, double temperature)
{
    // Reset phi (remove negative effects of irradiation)
    this->phi = 0;
    // Maybe change T1n?
    this->t1n *= 0.8;

    auto temp_tmp = this->system_temperature;
    this->set_system_temperature(temperature);
    this->run_for(t);
    this->set_system_temperature(temp_tmp);
}

double Simulation::find_optimal_freq(bool negative) const
{
    using boost::math::tools::brent_find_minima;

    // There's only a function to find the minimum, so we'll need to adjust our
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
    // TEMP_SS = steady-state temperature
    // K_TEMP = rate of exponential increase
    // If we're annealing, we shouldn't allow the temperature to change (assume
    // anneals occur at constant temperature)
    const double K_TEMP = 0.01;
    const double TEMP_SS = this->system_temperature + this->beam_current / 100;

    // Increase phi according to some exponential growth when the beam is on
    // Parameters are similar to those for temperature change
    const double K_PHI = this->beam_current / 1e7;
    const double PHI_SS = 0.001;

    // Calculate convenience constants
    const double A = this->t1e / this->t1n -
                     (this->c / 2) * (this->alpha + this->beta) - this->phi;
    const double B = (this->c / 2) * (this->alpha - this->beta);
    const double C = (this->alpha - this->beta) / 2;
    const double D = -1 - (this->alpha + this->beta) / 2;

    // Calculate rates
    const double PN_PRIME = (A * this->pn_raw + B * this->pe) / this->t1e;
    const double PE_PRIME =
        (C * this->pn_raw + D * this->pe + this->pe0) / this->t1e;

    // Update pn and pe using Euler's method
    this->pn_raw += PN_PRIME * t;
    this->pe += PE_PRIME * t;
    // Update temperature and phi
    this->set_temperature(this->temperature +
                          t * K_TEMP * (TEMP_SS - this->temperature));
    this->phi += t * K_PHI * (PHI_SS - this->phi);
    // Update C
    this->c += IRRADIATION_FACTOR * this->beam_current * t;

    // Update dose, along with the fit parameters M1 and M2, which depend on it.
    // Recall the exponential change of M1 and M2 described in the
    // documentation.
    const double DELTA_DOSE = (this->beam_current * 1e-9 / ELEM_CHARGE) * t;
    this->fit_params.m1 +=
        FIT_M1_COEFF * FIT_M1_RATE * DELTA_DOSE * exp(FIT_M1_RATE * this->dose);
    this->fit_params.m2 +=
        FIT_M2_COEFF * FIT_M2_RATE * DELTA_DOSE * exp(FIT_M2_RATE * this->dose);
    this->dose += (this->beam_current * 1e-9 / ELEM_CHARGE) * t;

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
    const double SCALE =
        this->fit_params.a / (sqrt(2 * PI) * this->fit_params.s);
    const double DIFF1 = freq - this->fit_params.m1;
    const double DIFF2 = freq - this->fit_params.m2;
    const double EXP1 =
        exp(-DIFF1 * DIFF1 / (2 * this->fit_params.s * this->fit_params.s));
    const double EXP2 =
        exp(-DIFF2 * DIFF2 / (2 * this->fit_params.s * this->fit_params.s));

    return std::make_pair(SCALE * EXP2, SCALE * EXP1);
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
