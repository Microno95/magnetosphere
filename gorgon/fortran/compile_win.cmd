call F2py -c streamtracer.f90 --f90flags='-fopenmp' -lgomp --fcompiler=gnu95 -m streamtracer
call F2py -c null_finder.f90 --f90flags='-fopenmp' -lgomp --fcompiler=gnu95 -m null_finder

::call F2py -c streamtracer.f90 --f90flags='/Qopenmp' -lgomp --fcompiler=intelem -m streamtracer
::call F2py -c null_finder.f90 --f90flags='/Qopenmp' -lgomp --fcompiler=intelem -m null_finder