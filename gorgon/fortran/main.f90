	include 'Streamtracer.f90'
    program main
    use streamtracer
    implicit none
    
    integer :: nx, ny, nz, nlines, i, j, k, il
    double precision, dimension(3) :: d
    double precision, dimension(:, :), allocatable :: x0
    double precision, dimension(:,:,:,:), allocatable :: vec
    double precision, dimension(:,:,:), allocatable :: xs
    integer, dimension(:,:,:), allocatable :: link
    double precision  :: x, y, z, r
    integer, dimension(:), allocatable ::  ROT, ns_out
    integer :: n_step
    
    nx = 200; ny = 100; nz = 100
    ns = 10000
    
    d = 10.*[2., 1., 1.]/[nx, ny, nz]
    xc = 0.5*d*[nx, ny, nz]
    
    n_step = 5
    
    nlines = nx*ny*nz/n_step**3
    
    allocate(vec(nx,ny,nz,3), link(nx,ny,nz), x0(nlines,3), xs(nlines, ns, 3))
    allocate(ROT(nlines), ns_out(nlines))
    
    il = 0
    do k=1,nz
        z = k*d(3) - xc(3)
        do j=1,ny
            y = j*d(2) - xc(2)
            do i=1,nx
                x = i*d(1) - xc(1)
                
                r = sqrt(x**2+y**2+z**2)
                
                vec(i,j,k,1) = (x/r)/r**3
                vec(i,j,k,2) = (y/r)/r**3
                vec(i,j,k,3) = (z/r-1)/r**3+1             
            end do
        end do
    end do
    
    il = 0
    do k=1,nz,n_step
        z = k*d(3) - xc(3)
        do j=1,ny,n_step
            y = j*d(2) - xc(2)
            do i=1,nx,n_step
                x = i*d(1) - xc(1)
                
                il = il+1
                x0(il,:) = [x, y, z]+xc
                
            end do
        end do
    end do
    
    ds = 0.1*d(1)
    inner_boundary=.true.
    call connectivity_array(x0, nlines, vec, nx, ny, nz, d, link)
    
    !print*,'Streamtracer', nlines
    call streamline_array(x0, nlines, vec, nx, ny, nz, d, 1, ns, xs, ROT, ns_out)
    
    end program main