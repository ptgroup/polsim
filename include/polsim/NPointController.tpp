namespace polsim
{
template <unsigned n_points>
NPointController<n_points>::NPointController(Pdp pdp, double step_size,
                                             bool seek_positive)
    : Controller(pdp, seek_positive), step_size(step_size)
{
}

template <unsigned n_points>
Data NPointController<n_points>::step()
{
    // Collect the specified number of data points.
    auto data = std::array<Data, n_points>();
    for (unsigned i = 0; i < data.size(); i++) {
        data[i] = pdp.take_data();
    }

    switch (make_decision(data)) {
    case Decision::SWITCH_DIRECTION:
        this->step_size *= STEP_SIZE_REDUCE;
        this->direction *= -1.0;
        // We shouldn't go below a certain step size.
        if (this->step_size < MIN_STEP_SIZE)
            this->step_size = MIN_STEP_SIZE;
        // fallthrough
    case Decision::KEEP_DIRECTION:
        this->pdp.set_freq(data.back().freq +
                           this->direction * this->step_size);
        break;
    case Decision::NO_MOTION:
        break;
    }

    return data.back();
}
} // namespace polsim
