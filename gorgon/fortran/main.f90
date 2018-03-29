    program main
    use connectivity_tracer
    implicit none
    
    integer :: nx, ny, nz, nlines, i, j, k, il
    double precision, dimension(3) :: d, xc
    double precision, dimension(:, :), allocatable :: x0
    double precision, dimension(:,:,:,:), allocatable :: vec
    integer, dimension(:,:,:), allocatable :: link
    double precision  :: x, y, z, r
    
    nx = 200; ny = 100; nz = 100
    
    d = 10.*[2., 1., 1.]/[nx, ny, nz]
    print*,d
    xc = 0.5*d*[nx, ny, nz]
    print*,xc
    
    nlines = nx*ny*nz/10**3
    print*,nlines
    
    allocate(vec(nx,ny,nz,3), link(nx,ny,nz), x0(nlines,3))
    
    il = 0
    do k=1,nz
        do j=1,ny
            do i=1,nx
                x = i*d(1) - xc(1)
                y = j*d(2) - xc(2)
                z = k*d(3) - xc(3)
                r = sqrt(x**2+y**2+z**2)
                
                vec(i,j,k,1) = (x/r)/r**3
                vec(i,j,k,2) = (y/r)/r**3
                vec(i,j,k,3) = (z/r-1)/r**3                
            end do
        end do
    end do
    
    il = 0
    do k=1,nz,10
        do j=1,ny,10
            do i=1,nx,10
                x = i*d(1) - xc(1)
                y = j*d(2) - xc(2)
                z = k*d(3) - xc(3)
                
                il = il+1
                x0(il,:) = [x, y, z]
                
            end do
        end do
    end do
    
    ns = 10000
    ds = 0.1*d(1)
    inner_boundary=.true.
    call connectivity_array(x0, nlines, vec, nx, ny, nz, d, xc, link)
    
    !k = floor(nz/2)
    !write(fmt,'("ES12.5,", I3, ) '
    !do j=1,ny
    !    write(100,fmt) (link(i,j,k), i=1,nx)
    !end do
    
    !ES12.5, nx(",", ES12.5)
    
    print*,minval(link)
    print*,maxval(link)
    
    end program main