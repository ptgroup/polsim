/**
 * @file pdp.hpp
 * @brief Simulation for the PDP data collection setup.
 *
 * @author Ian Johnson
 * @date June 29, 2017
 */
#ifndef PDP_HPP
#define PDP_HPP

#include "simulation.hpp"

#include <random>

namespace polsim
{
/// Time per sweep, in ms.
constexpr double MS_PER_SWEEP = 64;
/// Error in polarization per sweep.
constexpr double SWEEP_UNCERTAINTY = 0.01;

/**
 * @brief The PDP simulation.
 *
 * This class contains the state of the PDP simulation, which operates on an
 * underlying Simulation in a way analogous to the real system.
 */
class Pdp
{
	// The perfect controller knows everything
	friend class PerfectController;

	/// The underlying system.
	Simulation sim;
	/// The number of sweeps per data reading.
	unsigned n_sweeps;
	/// The internal random number generator.
	std::mt19937 rng;

public:
	/**
	 * @brief Constructs a new PDP simulator.
	 *
	 * @param sim The underlying simulation.
	 * @param n_sweeps The number of sweeps per data reading.
	 */
	Pdp(Simulation sim, unsigned n_sweeps = 250);

	/**
	 * @brief Gets a single data point by sweeping.
	 *
	 * As in the real system, this will take a certain number of data points
	 * with some uncertainty, and the returned data point will have a
	 * polarization value reflecting this. The rest of the data returned
	 * will
	 * correspond to the last sweep which was performed.
	 *
	 * @return The last data point taken, with the polarization obtained by
	 * sweeping.
	 */
	Data take_data();
	/**
	 * @brief Returns a reference to the underlying system.
	 *
	 * @return A reference to the underlying system.
	 */
	System &system_ref();
	/**
	 * @brief Sets the frequency of the underlying simulation.
	 *
	 * @param freq The new frequency.
	 */
	void set_freq(double freq);

private:
	/**
	 * @brief Performs a single sweep and returns its data.
	 */
	Data sweep();
};
}

#endif
