"""Build helper for the ni_dynamicclamp shared libraries.

Two extensions are built side by side:

* ``interface_c.so`` — legacy C backend. **FROZEN**: kept buildable for
  back-compat with ``ni_generic.py`` / ``lif_tester.py`` but no new
  features land here. See ``interface_c.c`` for details.
* ``libni.so`` — C++ backend mirroring ``interface.cpp`` (the same source
  Brian2's ``cpp_standalone`` codegen compiles inline). This is the
  shared library hot-loaded by ``ni_torch.py`` for the PyTorch path.

Most contributors should just run ``bash compile_cpp`` to rebuild
``libni.so``; this setup file exists so the build can be driven from
``pip``/``setuptools`` when packaging.
"""
from setuptools import Extension, setup
import os

# get the path to the current file
path = os.path.dirname(os.path.abspath(__file__))

setup(
    ext_modules=[
        Extension(
            name="interface_c.so",  # legacy C path — FROZEN
            include_dirs=[path],
            sources=["interface_c.c"],
            libraries=["nidaqmx"],
            extra_compile_args=["-g"],
        ),
        Extension(
            name="libni.so",  # C++ path used by ni_torch.py + Brian2 codegen
            include_dirs=[path],
            sources=["interface.cpp", "interface_torch_export.cpp"],
            libraries=["nidaqmx", "rt", "pthread"],
            extra_compile_args=["-O2", "-std=c++14", "-Wno-write-strings"],
            language="c++",
        ),
    ]
)