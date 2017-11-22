from distutils.core import setup

setup(name='Magnetosphere',
      version="0.0.1",
      description='Python subroutines for analysing Gorgon Output',
      author = "Lars Mejnertsen",
      author_email = "lars.mejnertsen10@imperial.ac.uk",
      packages=['boundary','models','vtk'],
      install_requires=['vtk', 'numpy', 'matplotlib', 'evtk'],
     )