// clang-format off

// clang-format on

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/types.hpp>
#include <cstdint>
#include <type_traits>

struct bubble_recorder
  {
    mutable std::map<double, uint64_t> weights;

    void operator()(const fwdpy11::singlepop_t& pop) const
      {
        double mkey;

        for (uint32_t mut = 0; mut < pop.mutations.size(); ++mut)
          {
	           mkey = pop.mutations[mut].pos;
             if (pop.mutations[mut].neutral)
              {
                if (pop.mcounts[mut] > 0)
                  {
                    if (weights.find(mkey) != weights.end())
                      {
                        weights[mkey] += pop.mcounts[mut];
                      }
                    else
		                  {
	                      weights[mkey] = pop.mcounts[mut];
                      }
                  }
              }
          }
      }
  };

namespace py = pybind11;

PYBIND11_PLUGIN(bubble_recorder)
{
  pybind11::module m("bubble_recorder",
		     "Record neutral mutant bubble sizes.");

  py::class_<bubble_recorder>(m, "BubbleRecorder")
    .def(py::init<>())
    .def_readonly("weights", &bubble_recorder::weights, "bubble weights")
    .def("__call__",
	 [](const bubble_recorder& r, const fwdpy11::singlepop_t& pop) -> void {
	   r(pop);
	 });

  return m.ptr();
}
