    module streamtracer
    use omp_lib
    implicit none
    
    ! Connectivity tracer!
    
    ! ROT: Reason of termination
    ! If  ROT = 0: Still running
    !     ROT = 1: Out of steps
    !     ROT = 2: Out of domain
    !     ROT = 3: In inner boundary
    !     ROT = -1: vmag = 0
    !     ROT = -2: NaN present
    ! Link: Connectivity in MS
    ! If  Link = 1: SW
    !     Link = 2: Closed
    !     Link = 3: North Open
    !     Link = 4: South Open
    !     Link = 5: SW Inc
    !     Link = 6: North Inc
    !     Link = 7: South Inc
    !     Link = 8: Inc Inc
    
    integer :: ns
    double precision :: ds, r_IB=1.
    double precision, dimension(3) :: xc
	logical :: inner_boundary=.false., write_threads=.false.
	
    contains
        
    subroutine streamline_array(x0, nlines, v, nx, ny, nz, d, dir, ns, xs, ROT, ns_out)
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns, nx, ny, nz, dir, nlines
    double precision, dimension(nlines, ns, 3), intent(out) :: xs
    integer, intent(out), dimension(nlines) :: ROT, ns_out
    double precision, dimension(3) :: x0_i
    double precision, dimension(ns, 3) :: xs_i
    integer :: ROT_i, ns_i
    integer :: i, j, k
	
	!if(write_threads) then
	!	!$omp parallel
	!	!$omp critical
	!	open(unit=500, file='streamtracer_threads.txt', access='append')
	!	write(500,*) omp_get_thread_num(), omp_get_max_threads()
	!	close(500)
	!	!$omp end critical
	!	!$omp end parallel
	!end if
    
	!$omp parallel do default(private) shared(xs,ROT,ns_out, x0) schedule(dynamic)
    do i=1,nlines
        x0_i = x0(i,:)
        call streamline(x0_i, v, nx, ny, nz, d, dir, ns, xs_i, ROT_i, ns_i)
        ROT(i) = ROT_i
        ns_out(i) = ns_i
        do j=1,ns_i
            xs(i,j,:) = xs_i(j,:)
        end do
        
    end do
	!$omp end parallel do
    
    end subroutine streamline_array

    subroutine streamline(x0, v, nx, ny, nz, d, dir, ns_in, xs, ROT, ns_out)
    implicit none
    double precision, dimension(3), intent(in) :: x0, d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, intent(in) :: ns_in, nx, ny, nz, dir
    double precision, dimension(ns_in, 3), intent(out) :: xs
    integer, intent(out) :: ROT, ns_out
    double precision, dimension(3) :: xi
    integer :: i
    
    ns = ns_in
    
    ROT = 0
    
    xs = 0
    xs(1, :) = x0
    xi = x0

    do i=2,ns
        
        ROT = check_bounds(xi, nx, ny, nz, d)
        if(ROT.ne.0) exit
        
        call RK4_update(xi, v, nx, ny, nz, d, dir)
        
        xs(i,:) = xi
        
    end do

    if (ROT.eq.0) ROT = 1
    ns_out = i-1
    
    end subroutine streamline
    
    subroutine connectivity_array(x0, nlines, v, nx, ny, nz, d, link)
    double precision, dimension(nlines,3), intent(in) :: x0
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    integer, dimension(nlines), intent(out) :: link
    integer, intent(in) :: nx, ny, nz, nlines
    integer :: ROT_f, ROT_r
    double precision, dimension(3) :: x0_i
    integer :: i
	
	!if(write_threads) then
	!	!$omp parallel
	!	!$omp critical
	!	open(unit=500, file='streamtracer_threads.txt', access='append')
	!	write(500,*) omp_get_thread_num(), omp_get_max_threads()
	!	close(500)
	!	!$omp end critical
	!	!$omp end parallel
	!end if
    
    !$omp parallel do default(firstprivate) shared(v, x0, link) private(x0_i) schedule(dynamic)
    do i=1,nlines
        x0_i = x0(i,:)
        ROT_f = streamline_end(x0_i, v, nx, ny, nz, d, 1)
        ROT_r = streamline_end(x0_i, v, nx, ny, nz, d, -1)
        link(i) = categorise_end_pts(ROT_f, ROT_r)
    end do
	!$omp end parallel do
    
    end subroutine connectivity_array
 
    integer function streamline_end(x0, v, nx, ny, nz, d, dir)
    implicit none
    double precision, dimension(3), intent(in) :: x0, d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    double precision, dimension(3) :: xi
    integer, intent(in) :: nx, ny, nz, dir
    integer :: ROT, i
    
    ROT = 0
    
    xi = x0

    do i=2,ns
        
        ROT = check_bounds(xi, nx, ny, nz, d)
        if(ROT.ne.0) exit
        
        call RK4_update(xi, v, nx, ny, nz, d, dir)
        
    end do
 
    if (ROT.eq.0) ROT = 1
    
    streamline_end = ROT
 
    end function streamline_end
    
    integer function categorise_end_pts(f, r)
    integer, intent(in) :: f, r
    integer :: link
    
    if(r==2 .and. f==2) then
		link = 1	! Solar Wind
	elseif(r==3 .and. f==3) then
		link = 2	! Closed
	elseif(r==3 .and. f==2) then
		link = 3	! North-Open
	elseif(r==2 .and. f==3) then
		link = 4	! South-Open
	elseif(r==2 .or. f==2) then
		link = 5	! SW-Inc
	elseif(r==3) then
		link = 6	! North-Inc
	elseif(f==3) then
		link = 7	! South-Inc
	else
		link = 8	! Inc-Inc
    end if  
            
    categorise_end_pts = link
    
    end function categorise_end_pts
    
    subroutine RK4_update(xi, v, nx, ny, nz, d, dir)
    double precision, dimension(3), intent(inout) :: xi
    double precision, dimension(3), intent(in) :: d
    double precision, dimension(nx,ny,nz,3), intent(in) :: v
    double precision, dimension(3) :: xu
    integer, intent(in) :: nx, ny, nz, dir
    double precision, dimension(3) :: k1, k2, k3, k4
    integer :: i

    !--- RK4 K parameters ---------------------------------------------------------------------
    call stream_function(xi, v, nx, ny, nz, d, dir, k1)
        
    xu = xi+0.5*k1
    call stream_function(xu, v, nx, ny, nz, d, dir, k2)
        
    xu = xi+0.5*k2
    call stream_function(xu, v, nx, ny, nz, d, dir, k3)
        
    xu = xi+k3
    call stream_function(xu, v, nx, ny, nz, d, dir, k4)

    !--- Step ---------------------------------------------------------------------------------
        
    xi = xi + (k1 + 2*k2 + 2*k3 + k4)/6.
    
    end subroutine RK4_update
    
    double precision function check_bounds(xi, nx, ny, nz, d)
    double precision, intent(in), dimension(3) ::xi, d
    integer, intent(in) :: nx, ny, nz
    double precision :: ri
    
    check_bounds = 0
    
    if ( isnan(xi(1)).or.isnan(xi(2)).or.isnan(xi(3)) ) then
        check_bounds = -2
        return
    end if
    
    if(xi(1).lt.0.or.xi(1).gt.d(1)*nx) then
        check_bounds = 2
        return
    end if
    
    if(xi(2).lt.0.or.xi(2).gt.d(2)*ny) then
        check_bounds = 2
        return
    end if
    
    if(xi(3).lt.0.or.xi(3).gt.d(3)*nz) then
        check_bounds = 2
        return
    end if
    
    
    if(inner_boundary) then
        ri = sqrt((xi(1)-xc(1))**2+(xi(2)-xc(2))**2+(xi(3)-xc(3))**2)
        if(ri.le.r_IB) then
            check_bounds = 3
            return
        end if
    end if
    
    end function check_bounds

    subroutine stream_function(xI, v, nx, ny, nz, d, dir, f)
    implicit none
    double precision, dimension(3), intent(in) :: xI
    double precision, dimension(nx, ny, nz, 3), intent(in) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(3), intent(in) :: d
    integer, intent(in) :: dir
    double precision, dimension(3), intent(out) :: f
    double precision :: vmag
    double precision, dimension(3) :: vI
    
    call interpolate(xI, v, nx, ny, nz, d, vI)
    
	vmag  = sqrt(vI(1)**2+vI(2)**2+vI(3)**2)
    f = dir*vI/vmag*ds
	
    end subroutine stream_function
    
    subroutine interpolate(xI, v, nx, ny, nz, d, vI)
    double precision, intent(in), dimension(3) :: xI
    double precision, intent(in), dimension(nx,ny,nz,3) :: v
    integer, intent(in) :: nx, ny, nz
    double precision, intent(in), dimension(3) :: d
    double precision, intent(out), dimension(3) :: vI
    double precision, dimension(3) :: distI
    integer, dimension(3) :: i0, i1
    double precision, dimension(2,2,2) :: cell
    
    i0 = floor(xI/d)
    i0(1) = min(max(1,i0(1)), nx-1)
    i0(2) = min(max(1,i0(2)), ny-1)
    i0(3) = min(max(1,i0(3)), nz-1)
    
    i1 = i0+1
    
    distI = xI/d-i0
    
    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 1)
    call interp_trilinear(distI, cell, vI(1))
    
    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 2)
    call interp_trilinear(distI, cell, vI(2))
    
    cell = v(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), 3)
    call interp_trilinear(distI, cell, vI(3))
    
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

    end module streamtracer