/**
 * @file pdp.hpp
 * @brief Simulation for the PDP data collection setup.
 *
 * @author Ian Johnson
 * @date June 29, 2017
 */
#ifndef _PDP_HPP
#define _PDP_HPP

#include "simulation.hpp"

#include <random>

/// Time per sweep, in ms.
constexpr double MS_PER_SWEEP = 64;
/// Error in polarization per sweep.
constexpr double SWEEP_UNCERTAINTY = 0.04;

/**
 * @brief The PDP simulation.
 *
 * This class contains the state of the PDP simulation, which operates on an
 * underlying Simulation in a way analogous to the real system.
 */
class Pdp
{
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
    Pdp(Simulation &&sim, unsigned n_sweeps = 250);

    /**
     * @brief Gets a single data point by sweeping.
     *
     * As in the real system, this will take a certain number of data points
     * with some uncertainty, and the returned data point will have a
     * polarization value reflecting this. The rest of the data returned will
     * correspond to the last sweep which was performed.
     *
     * @return The last data point taken, with the polarization obtained by
     * sweeping.
     */
    Data take_data();

  private:
    /**
     * @brief Performs a single sweep and returns its data.
     */
    Data sweep();
};

#endif
