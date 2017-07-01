#include "pdp.hpp"

Pdp::Pdp(Simulation sim, unsigned n_sweeps) : sim(sim), n_sweeps(n_sweeps)
{
    this->rng.seed(std::random_device{}());
}

Data Pdp::take_data()
{
    // The average of the polarization readings obtained by sweeping
    auto pn = 0.0;
    for (auto i = 0U; i < this->n_sweeps; i++)
        pn += this->sweep().pn;
    pn /= this->n_sweeps;

    // Get a data point and change its polarization
    auto data = this->sim.take_data();
    data.pn = pn;
    return data;
}

System &Pdp::system_ref() { return this->sim.system_ref(); }

void Pdp::set_freq(double freq) { this->sim.set_freq(freq); }

Data Pdp::sweep()
{
    this->sim.run_for(MS_PER_SWEEP / 1000);
    auto sweep_data = this->sim.take_data();

    // Fuzz polarization
    auto dist = std::uniform_real_distribution<double>{-1, 1};
    sweep_data.pn += SWEEP_UNCERTAINTY * dist(rng);
    return sweep_data;
}
