// clang-format off
<% 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'] = [ fp11.get_includes(), fp11.get_fwdpp_includes() ]
#On OS X using clang, there is more work to do.  Using gcc on OS X
#gets rid of these requirements. The specifics sadly depend on how
#you initially built fwdpy11, and what is below assumes you used
#the provided setup.py + OS X + clang:
#cfg['compiler_args'].extend(['-stdlib=libc++','-mmacosx-version-min=10.7'])
#cfg['linker_args']=['-stdlib=libc++','-mmacosx-version-min=10.7']
setup_pybind11(cfg)  
#An alternative to the above is to add the first line to CPPFLAGS
#and the second to LDFLAGS when compiling a plugin on OS X using clang.
%>
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
