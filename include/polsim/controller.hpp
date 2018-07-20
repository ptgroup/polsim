/**
 * @file controller.hpp
 * @brief Simulations for various motor controller algorithms.
 *
 * @author Ian Johnson
 * @date June 30, 2017
 */
#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include "polsim/pdp.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <list>
#include <numeric>
#include <utility>
#include <vector>

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
    virtual ~Controller() = default;

    /**
     * @brief Performs a single step of the algorithm.
     *
     * Usually, what this will look like is that the controller will take
     * several data points in order to make a decision concerning how to move to
     * motor, and will then move the motor, concluding the step. This behavior
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
 * @brief A controller that uses a fixed number of data points on each step.
 *
 * @tparam n_points The number of data points to take per step.
 */
template <unsigned n_points>
class NPointController : public Controller
{
    /// The minimum step size.
    constexpr static double MIN_STEP_SIZE = 0.001;
    /// The fraction by which to decrease step size on direction change.
    constexpr static double STEP_SIZE_REDUCE = 0.8;
    /// The initial step size (GHz).
    double initial_step_size;
    /// The frequency step size (GHz).
    double step_size;
    /// The direction (the sign of the step size).
    double direction = 1.0;

protected:
    enum class Decision { NO_MOTION, KEEP_DIRECTION, SWITCH_DIRECTION };

    /**
     * @brief Given the data points for the current step, return what to do
     * next.
     */
    virtual Decision make_decision(const std::array<Data, n_points> &data) = 0;

public:
    /**
     * @brief Constructs a new controller.
     *
     * @param pdp The underlying PDP simulation.
     * @param step_size The initial step size to use (in GHz).
     * @param seek_positive Whether to seek positive polarization.
     */
    NPointController(Pdp pdp, double step_size = 0.05,
                     bool seek_positive = true);

    /**
     * @brief Resets the step size.
     */
    virtual void reseek();
    Data step() override;
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
 * 3).
 */
template <unsigned n_points>
class StandardController : public NPointController<n_points>
{
    static_assert(n_points >= 3, "Must use at least 3 points per step");

    /// The time (in seconds) at which to switch algorithms.
    constexpr static double ALGO_SWITCH_TIME = 600.0;
    /// The "good" rate ratio threshold.
    constexpr static double GOOD_RATIO = 0.8;
    /// The number of bad polarization values to accept before reseeking.
    constexpr static unsigned RESEEK_AFTER = 4;
    /// The polarization spread to be considered "steady".
    constexpr static double STEADY_POL_SPREAD = 0.002;

    /// The "countdown timer" (in seconds) to switch algorithms to polarization.
    double algo_switch_countdown = ALGO_SWITCH_TIME;
    /// The number of bad polarization values collected so far.
    unsigned bad_pols = 0;
    /// The last k value calculated.
    double last_k = 0.0;
    /// The last polarization value collected.
    double last_pol = 0.0;
    /// The time of the last data point.
    double last_time = 0.0;
    /// Whether to seek using polarization values.
    bool seek_pol = false;

protected:
    typename NPointController<n_points>::Decision
    make_decision(const std::array<Data, n_points> &data) override;

    /// Makes the decision based on rate values.
    typename NPointController<n_points>::Decision
    make_decision_rate(const std::array<Data, n_points> &data);

    /// Makes the decision based on polarization values.
    typename NPointController<n_points>::Decision
    make_decision_polarization(const std::array<Data, n_points> &data);

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

    void reseek() override;
    Data step() override;
};

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
} // namespace polsim

#include "NPointController.tpp"
#include "StandardController.tpp"

#endif
