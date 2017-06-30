/**
 * @file simulation.hpp
 * @brief The simulation for the solid polarized target.
 *
 * @author Ian Johnson
 * @date June 29, 2017
 */
#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#include <iostream>
#include <random>

/// Elementary charge (in C).
constexpr double ELEM_CHARGE = 1.602176662e-19;
/// Pi.
constexpr double PI = 3.1415926;

/**
 * @brief The rate at which the fit parameter M1 changes with respect to dose.
 *
 * The M1 parameter changes exponentially with dose:
 * `M1 = C1 + ::FIT_M1_COEFF * exp(::FIT_M1_RATE * dose)`
 * where `C1` is some constant.
 */
constexpr double FIT_M1_RATE = -0.38e-15;
/// See the documentation for ::FIT_M1_RATE.
constexpr double FIT_M1_COEFF = 0.045;
/// The same as ::FIT_M1_RATE, but for M2.
constexpr double FIT_M2_RATE = -3.8e-15;
/// The same as ::FIT_M1_COEFF, but for M2.
constexpr double FIT_M2_COEFF = -0.065;

/**
 * @brief Thermal randomness fraction.
 *
 * Due to temperature fluctuations (not actually simulated), the polarization
 * will fluctuate slightly. This parameter gives the fraction bounding the
 * fluctuations; that is, it will be multiplied by the current polarization to
 * give the maximum deviation.
 */
constexpr double THERMAL_RANDOMNESS = 0.02;
/**
 * @brief Base randomness.
 *
 * This is an absolute polarization amount; the polarization will experience
 * fluctuations of this sort no matter what, and the amount of fluctuation is
 * NOT dependent on the current polarization (unlike ::THERMAL_RANDOMNESS).
 */
constexpr double BASE_RANDOMNESS = 0.002;

/**
 * @brief Irradiation factor (in 1 / (nA * s))
 *
 * Multiplied by the beam current, this gives the rate of change of C over time
 * (C will drift when the beam is on).
 */
constexpr double IRRADIATION_FACTOR = 1e-10;

/**
 * @brief The fit parameters which are used to calculate alpha and beta.
 *
 * Alpha and beta are modelled as Gaussian distributions with respect to
 * frequency. The default values are taken from a fit performed on data from
 * the April 2016 cooldown.
 */
struct FitParameters {
    /// The height of the distributions.
    double a = 0.545266;
    /// The standard deviation.
    double s = 0.088415;
    /// The mean of the beta distribution.
    double m1 = 140.241;
    /// The mean of the alpha distribution.
    double m2 = 140.533;
};

/**
 * @brief A single observable data point.
 *
 * Each data point contains all observable data at a particular time.
 */
struct Data {
    double t, pn, pe, freq, c, temperature, dose;

    /**
     * @brief Outputs the data in standard CSV format.
     *
     * This will be `time, pn, freq`.
     */
    friend std::ostream &operator<<(std::ostream &out, const Data &data);
};

/**
 * @brief The simulation for spin 1/2.
 */
class Simulation
{
    /// The fit parameters currently in use.
    FitParameters fit_params;

    double t1n = 25 * 60;
    double t1e = 0.03;

    /// Current time (in seconds).
    double t = 0;
    /// Frequency (GHz).
    double freq;
    /// Temperature (K).
    double temperature;
    /// System temperature (K).
    double system_temperature = 1;

    double alpha, beta;
    double c = 0.000136073;
    double pe0;
    double phi = 0;

    /// The "raw polarization" (no random fluctuations).
    double pn_raw = 0;
    double pn = 0;
    double pe = -1;

    /// Dose (electrons per cm^3).
    double dose = 0;
    /// Beam current (nA).
    double beam_current = 0;

    /// The internal random number generator.
    std::mt19937 rng;

  public:
    /**
     * @brief Constructs a Simulation.
     *
     * @param freq The initial frequency.
     */
    Simulation(double freq);

    /**
     * @brief Returns observable data readings.
     */
    Data take_data();
    /**
     * @brief Sets the frequency.
     *
     * @param freq The new frequency.
     */
    void set_freq(double freq);
    /**
     * @brief Sets the system temperature.
     *
     * The system temperature is distinct from the internal temperature used in
     * the simulation. The system temperature may be thought of as the
     * temperature of the fridge, while the internal temperature is the
     * temperature of the target itself.
     *
     * @param temperature The new system temperature.
     */
    void set_system_temperature(double temperature);
    /**
     * @brief Turns the beam on.
     *
     * @param current The beam current, in nA.
     */
    void beam_on(double current);
    /**
     * @brief Turns the beam off.
     */
    void beam_off();
    /**
     * @brief Runs the simulation for the specified time.
     *
     * Note that if the given time is not a multiple of the time step, the
     * simulation will actually run for longer than the specified time. Large
     * time steps will lead to more inaccuracy in the simulated values.
     *
     * @param t The time (in seconds) to run the simulation.
     * @param step The time step to use (as in Euler's method).
     */
    void run_for(double t, double step = 0.001);
    /**
     * @brief Performs an anneal.
     *
     * The mechanics of the anneal are not finalized, and may be completely
     * bogus.
     *
     * @param t The time to anneal (in seconds).
     * @param temperature The temperature at which to anneal (in K).
     */
    void anneal(double t, double temperature);

  private:
    /**
     * @brief Sets the internal simulation temperature.
     *
     * This ensures that any dependent parameters are set correctly.
     */
    void set_temperature(double temperature);
    /**
     * @brief Makes a single time step.
     *
     * This performs a single iteration of Euler's method using the solid effect
     * equations, where the step size is the specified time.
     *
     * @param t The time step to use.
     */
    void time_step(double t);
    /**
     * @brief Calculates the parameters alpha and beta.
     */
    void calc_transition_rates();
    /**
     * @brief Returns the polarization after taking into account thermal
     * fluctuations.
     */
    double pn_noisy();
};

#endif
