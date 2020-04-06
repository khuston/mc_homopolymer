#include <pybind11/pybind11.h>
#include <mc_chain.hpp>

namespace py = pybind11;

using namespace Polymers;

PYBIND11_MODULE(PyPolymers, module) {
    module.doc() = "Polymers module";
    py::class_<MonteCarloChain>(module, "MonteCarloChain")
        .def(py::init<const int&, const double&, const double&>())
        .def("run", &MonteCarloChain::run)
        .def("__repr__",
            [](const MonteCarloChain &a) {
                return "<Polymers.MonteCarloChain>";
            }
        );
}