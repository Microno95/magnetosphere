    module connectivity_tracer
    use omp_lib
    implicit none
    
    ! Connectivity tracer!
    
    ! ROT: Reason of termination
    ! If  ROT = 0: Still running
    !     ROT = 1: Out of steps
    !     ROT = 2: Out of domain
    !     ROT = 3: In inner boundary: North
    !     ROT = 4: In inner boundary: South
    !     ROT = -1: vmag = 0
    !     ROT = -2: NaN present
    ! Link: Connectivity in MS
    ! If  Link = 0: SW
    ! If  Link = 1: Closed
    ! If  Link = 2: North Open
    ! If  Link = 3: South Open
    
    integer :: ns
    double precision :: ds
	
    contains
    
    subroutine streamline_array(x0, nlines, v, nx, ny, nz, d, xc, link)
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d, xc
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, dimension(nlines), intent(out) :: link
    integer, intent(in) :: nx, ny, nz, nlines
    integer :: ROT_f, ROT_r
    integer :: i
    
    !$omp parallel do default(firstprivate) shared(v, x0, link) schedule(dynamic)
    do i=1,nlines
        ROT_f = streamline(x0(i,:), v, nx, ny, nz, d, xc, 1)
        ROT_r = streamline(x0(i,:), v, nx, ny, nz, d, xc, -1)
        link(i) = categorise_end_pts(ROT_f, ROT_r)
    end do
	!$omp end parallel do
    
    end subroutine streamline_array
    
    integer function categorise_end_pts(f, r)
    integer, intent(in) :: f, r
    integer :: link
    
    if(f==2 .and. r==2) then
        link = 1
    elseif( (f==3 .or. f==4) .and. (r==3 .or. r==4) ) then
        link = 2
    elseif( (f==2 .and. r==3) .or. (f==3 .and. r==2) ) then
        link = 3
    elseif( (f==2 .and. r==4) .or. (f==4 .and. r==2) ) then
        link = 4
    elseif( f==2 .or. r==2) then
        link = 5
    elseif( f==3 .or. r==3) then
        link = 6
    elseif( f==4 .or. r==4) then
        link = 7
    else
        link = -1
    end if  
            
    categorise_end_pts = link
    
    end function categorise_end_pts

    integer function streamline(x0, v, nx, ny, nz, d, xc, dir)
    implicit none
    double precision, dimension(3), intent(in) :: x0, d, xc
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    double precision, dimension(3) :: xi
    integer, intent(in) :: nx, ny, nz, dir
    integer :: ROT
    double precision, dimension(3) :: k1, k2, k3, k4
    integer :: i
    
    ROT = 0
    
    xi = x0

    do i=0,ns-1

        !--- RK4 K parameters ---------------------------------------------------------------------
        call stream_function(xi,        v, nx, ny, nz, d, xc, dir, k1, ROT)
        k1 = k1*ds

        if(ROT.ne.0) exit
        
        call stream_function(xi+0.5*k1, v, nx, ny, nz, d, xc, dir, k2, ROT)
        k2 = k2*ds

        if(ROT.ne.0) exit
        
        call stream_function(xi+0.5*k2, v, nx, ny, nz, d, xc, dir, k3, ROT)
        k3 = k3*ds

        if(ROT.ne.0) exit
        
        call stream_function(xi+k3,     v, nx, ny, nz, d, xc, dir, k4, ROT)
        k4 = k4*ds

        if(ROT.ne.0) exit

        !--- Step ---------------------------------------------------------------------------------

        xi = xi + (k1 + 2*k2 + 2*k3 + k4)/6.

        if ( isnan(xi(1)).or.isnan(xi(2)).or.isnan(xi(3)) ) then
            ROT = -2
            exit
        end if
        
    end do

    if (ROT.eq.0) ROT = 1
    
    streamline = ROT

    end function streamline

    subroutine stream_function(xI, v, nx, ny, nz, d, xc, dir, f, ROT)
    implicit none
    double precision, dimension(3), intent(in) :: xI, xc
    double precision, dimension(nx, ny, nz, 3), intent(in) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(3), intent(in) :: d
    integer, intent(in) :: dir
    double precision, dimension(3), intent(out) :: f
    integer, intent(out) :: ROT
    double precision :: ri, vmag
    double precision, dimension(3) :: vI
    
    ri = sqrt(xI(1)**2 + xI(2)**2 + xI(3)**2)
    if(ri.lt.1.) then
        if(xI(3).ge.0) ROT = 3
        if(xI(3).lt.0) ROT = 4
        return
    end if
    
    call interpolate(xI, v, nx, ny, nz, d, xc, vI, ROT)
    
    if(ROT.ne.0) return
	
	vmag  = sqrt(vI(1)**2+vI(2)**2+vI(3)**2)
	if(vmag.eq.0.) then
		ROT = -1
		f = [0., 0., 0.]
	else
		f = dir*vI/vmag
    end if
	
    end subroutine stream_function
    
    subroutine interpolate(xI, v, nx, ny, nz, d, xc, vI, ROT)
    double precision, intent(in), dimension(3) :: xI, xc
    double precision, intent(in), dimension(nx,ny,nz,3) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, intent(in), dimension(3) :: d
    double precision, intent(out), dimension(3) :: vI
    integer, intent(out) :: ROT
    double precision, dimension(3) :: distI
    integer, dimension(3) :: i0, i1, n
    
    n = (/nx, ny, nz/)
    i0 = floor((xI + xc)/d)
    i1 = i0+1
    
    distI = (xI + xc)/d-i0
    
    !print*, 'xi: ', xi
    !print*, 'ixi: ', d*i0
    !print*, 'd: ', d
    !print*, 'i: ', i0
    !print*, 'n: ', n
    !print*, 'dist: ', disti
    
    if(i0(1).lt.1.or.i1(1).gt.nx-1.or. &
       i0(2).lt.1.or.i1(2).gt.ny-1.or. &
       i0(3).lt.1.or.i1(3).gt.nz-1) then
        ROT = 2
        vI = (/0., 0., 0./)
        return
    end if
    
    call interp_trilinear(distI, v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 1), vI(1))
    call interp_trilinear(distI, v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 2), vI(2))
    call interp_trilinear(distI, v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 3), vI(3))
    
    !print*, 'vI: ', vI
    
    end subroutine interpolate
    

    !--- Trilinear interpolation function -------------------------------------------------------------

    subroutine interp_trilinear(xd, f, fI)
    implicit none
    double precision, intent(in), dimension(3) :: xd
    double precision, intent(in), dimension(0:1,0:1,0:1) :: f
    double precision, intent(out) :: fI
    double precision, dimension(0:1,0:1) :: c
    double precision, dimension(3) :: m_xd
    double precision :: c0, c1
    
    m_xd = 1-xd

    !--- Interpolate over x -----------------------------------------------------------------------

    c(0,0) = f(0,0,0)*m_xd(1) + f(1,0,0)*xd(1)
    c(1,0) = f(0,1,0)*m_xd(1) + f(1,1,0)*xd(1)
    c(0,1) = f(0,0,1)*m_xd(1) + f(1,0,1)*xd(1)
    c(1,1) = f(0,1,1)*m_xd(1) + f(1,1,1)*xd(1)

    !--- Interpolate over y -----------------------------------------------------------------------

    c0 = c(0,0)*m_xd(2) + c(1,0)*xd(2)
    c1 = c(0,1)*m_xd(2) + c(1,1)*xd(2)

    !--- Interpolate over z -----------------------------------------------------------------------

    fI = c0*m_xd(3) + c1*xd(3)

    end subroutine interp_trilinear
    
    !--- RK4 streamline integrator and interpolation --------------------------------------------------

    end module connectivity_tracer