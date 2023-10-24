from setuptools import Extension, setup
import os
#get the path to the current file
path = os.path.dirname(os.path.abspath(__file__))
setup(
    ext_modules=[
        Extension(
            name="interface_c.so",  
            include_dirs = [path],
            sources=["interface_c.c"], # all sources are compiled into a single binary file
            libraries=["nidaqmx"],
            extra_compile_args=["-g"]
        ),
    ])