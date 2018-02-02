call F2py -c connectivity.f90 --f90flags='/Qopenmp' -m connectivity_tracer
call F2py -c Streamtracer.f90 --f90flags='/Qopenmp' -m streamtracer
call F2py -c null_finder.f90 --f90flags='/Qopenmp' -m null_finder