call F2py -c streamtracer.f90 --f90flags='/Qopenmp' --fcompiler=intelvem -m streamtracer
call F2py -c null_finder.f90 --f90flags='/Qopenmp' --fcompiler=intelvem -m null_finder