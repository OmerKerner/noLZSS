#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
#include <string>
#include <string_view>
#include <stdexcept>
#include "factorizer.hpp"
#include "version.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_noLZSS, m) {
    m.doc() = "Non-overlapping Lempel-Ziv-Storer-Szymanski factorization";

    py::class_<noLZSS::Factor>(m, "Factor")
        .def_readonly("start", &noLZSS::Factor::start)
        .def_readonly("length", &noLZSS::Factor::length);

    m.def("factorize", [](py::buffer b) {
        // Accept any bytes-like 1-byte-per-item contiguous buffer (e.g. bytes, bytearray, memoryview)
        py::buffer_info info = b.request();
        if (info.itemsize != 1) {
            throw std::invalid_argument("factorize: buffer must be a bytes-like object with itemsize==1");
        }
        if (info.ndim != 1) {
            throw std::invalid_argument("factorize: buffer must be a 1-dimensional bytes-like object");
        }

        const char* data = static_cast<const char*>(info.ptr);
        std::string_view sv(data, static_cast<size_t>(info.size));

        // Release GIL while doing heavy C++ work
        py::gil_scoped_release release;
        auto factors = noLZSS::factorize(sv);
        py::gil_scoped_acquire acquire;

        py::list out;
        for (auto &f : factors) out.append(py::make_tuple(f.start, f.length));
        return out;
    }, py::arg("data"), R"doc(Return list of (start,length) factors. Accepts bytes-like 1-D buffer.)doc");

    m.attr("__version__") = std::to_string(noLZSS::VERSION_MAJOR) + "." + std::to_string(noLZSS::VERSION_MINOR) + "." + std::to_string(noLZSS::VERSION_PATCH);
}
