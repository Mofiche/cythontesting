from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize(["cythontesting/water_saturation.pyx",
                             "cythontesting/say_hello.pyx"]))