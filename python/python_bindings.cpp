#include <pybind11/pybind11.h>
#include <mc_chain.h>

namespace py = pybind11;

using namespace Polymers;

PYBIND11_MODULE(Polymers, module) {
    module.doc() = "Polymers module";
    py::class_<MonteCarloChain>(module, "MCChain")
        .def(py::init<const int&, const double&, const double&>())
        .def("run", &MonteCarloChain::run)
        .def("__repr__",
            [](const MonteCarloChain &a) {
                return "<Polymers.MonteCarloChain>";
            }
        );        
}