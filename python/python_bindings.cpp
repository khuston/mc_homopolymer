#include <pybind11/pybind11.h>
#include <mc_chain.hpp>
#include <gsl/gsl_rng.h>

namespace py = pybind11;

using namespace Polymers;

PYBIND11_MODULE(Polymers, module) {
    module.doc() = "Polymers module";
    py::class_<MonteCarloChain>(module, "MonteCarloChain")
        .def(py::init<const int&, const double&, const double&, const unsigned long int&>(), py::arg("count"), py::arg("epsilon"), py::arg("sigma"), py::arg("seed") = gsl_rng_default_seed)
        .def("run", &MonteCarloChain::run, py::arg("nsteps"), py::arg("tethered"))
        .def("__repr__",
            [](const MonteCarloChain &a) {
                return "<Polymers.MonteCarloChain>";
            }
        );
}