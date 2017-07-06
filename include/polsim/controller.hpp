/**
 * @file controller.hpp
 * @brief Simulations for various motor controller algorithms.
 *
 * @author Ian Johnson
 * @date June 30, 2017
 */
#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include "pdp.hpp"

namespace polsim
{
/**
 * @brief The basic (abstract) representation of a motor controller.
 *
 * All motor controllers must derive from this class, which provides only the
 * most basic definitions needed to create a consistent interface.
 */
class Controller
{
protected:
	/// The underlying PDP simulation, from which all data is taken.
	Pdp pdp;
	/// Whether to seek positive polarization.
	bool seek_positive;

public:
	/**
	 * @brief Constructs a controller.
	 *
	 * @param pdp The underlying PDP simulation.
	 * @param seek_positive Whether to seek positive polarization.
	 */
	Controller(Pdp pdp, bool seek_positive = true);

	/**
	 * @brief Performs a single step of the algorithm.
	 *
	 * Usually, what this will look like is that the controller will take
	 * several data points in order to make a decision concerning how to
	 * move to
	 * motor, and will then move the motor, concluding the step. This
	 * behavior
	 * may differ, though, depending on the controller.
	 *
	 * @return The last data point collected by the underlying PDP.
	 */
	virtual Data step() = 0;
	/**
	 * @brief Returns a reference to the underlying system.
	 *
	 * @return A reference to the underlying system.
	 */
	System &system_ref();
};

/**
 * @brief The "standard" motor controller.
 *
 * This motor controller uses a simple "rate comparison" algorithm; that is, it
 * simply collects a bunch of data points, averages the rate of polarization
 * change over time, and compares that rate to the one previously collected. If
 * the rate goes down (for positive seek), then it switches direction, with a
 * small decrease in step size. If the rate goes up, it keeps going in the same
 * direction.
 *
 * @tparam n_points The number of data points to take per step (must be at least
 * 2).
 */
template <unsigned n_points>
class StandardController : public Controller
{
	static_assert(n_points >= 2, "Must use at least 2 points per step");

	/// The fraction by which to decrease step size on direction change.
	constexpr static double STEP_SIZE_REDUCE = 0.8;

	/// The last rate collected by this controller.
	double last_rate = 0;
	/// The frequency step size (GHz).
	double step_size;

public:
	/**
	 * @brief Constructs a new controller.
	 *
	 * @param pdp The underlying PDP simulation.
	 * @param step_size The initial step size to use (in GHz).
	 * @param seek_positive Whether to seek positive polarization.
	 */
	StandardController(Pdp pdp, double step_size = 0.05,
	                   bool seek_positive = true);

	Data step() override;
};

#include "StandardController.tpp"

/**
 * @brief A perfect motor controller.
 *
 * This motor controller uses knowledge of the optimal frequency at every point
 * to conduct a perfect seek. Note that this is not possible in the real world
 * without making the assumption that every system will behave exactly like this
 * simulation (and it also uses some parameters that would not be measurable in
 * any case).
 */
class PerfectController : public Controller
{
public:
	PerfectController(Pdp pdp, bool seek_positive = true);

	Data step() override;
};

/**
 * @brief An experimental motor controller.
 */
class ExperimentalController : public Controller
{
public:
	ExperimentalController(Pdp pdp, bool seek_positive = true);

	Data step() override;
};
}

#endif
