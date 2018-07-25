# polsim

This is a simulation for the solid polarized target experiment.

## Building

To build the project, you will need the following dependencies installed:

- A fairly modern C++ compiler (with C++14 support)
- The Boost C++ math library
- The Meson build system and the Ninja build tool
- Doxygen (optional, for building documentation)

On Debian/Ubuntu systems, this can be done by running the following command
as root:

```sh
apt install g++ libboost-math-dev meson doxygen
```

Then, you can run the following command to build everything:

```sh
mkdir build && meson build && ninja -C build
```

The `polsim` executable will be located in `build/src/polsim`, and the
documentation (if Doxygen was installed) can be viewed by opening the
`build/docs/html/index.html` file in a web browser.

## Modifying the simulated behavior

Right now, there is no support for running external "scripts" for simulating
certain scenarios, as previous simulations have done. Instead, you should
modify `src/main.cpp` directly with what you want. You can refer to the
generated documentation for the available methods.
