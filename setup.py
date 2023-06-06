from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['spkmeans.c', 'spkmeansmodule.c', 'utils.c'])

setup(name='mykmeanssp',
      version='1.0',
      description='Python wrapper for custom C extension',
      ext_modules=[module])
