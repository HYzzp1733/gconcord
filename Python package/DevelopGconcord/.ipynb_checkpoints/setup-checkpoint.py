from distutils.core import setup

setup (name = "gconcord",
       version = '1.0.0',
       description = "Python package of graphical CONCORD.",
       author = "Zhipu Zhou",
       packages = ['gconcord'],
       package_dir = {'gconcord':'gconcord'},
       package_data = {'gconcord':['sharedlib.so']},
      )