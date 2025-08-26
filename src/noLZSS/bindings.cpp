#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
#include "factorizer.hpp"
#include "version.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_noLZSS, m) {
    m.doc() = "Non-overlapping Lempel–Ziv–Storer–Szymanski factorization";

    py::class_<noLZSS::Factor>(m, "Factor")
        .def_readonly("start", &noLZSS::Factor::start)
        .def_readonly("length", &noLZSS::Factor::length);

    m.def("factorize", [](py::buffer b) {
        py::buffer_info info = b.request();
        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, info.size);
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize(sv);
        py::gil_scoped_acquire acquire;
        py::list out; out.reserve(factors.size());
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length));
        return out;
    }, py::arg("data"), R"doc(Return list of (start,length) factors.)doc");

    m.attr("__version__") = py::format("{}.{}.{}", noLZSS::VERSION_MAJOR, noLZSS::VERSION_MINOR, noLZSS::VERSION_PATCH);
}
