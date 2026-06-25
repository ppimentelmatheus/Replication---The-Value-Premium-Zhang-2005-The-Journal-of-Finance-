subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!--------------------------------------------------------------
!                 MEX file for SSIEFCN routine
!--------------------------------------------------------------
implicit none

integer plhs(*), prhs(*)
integer nlhs, nrhs
integer mxGetM, mxGetN, mxGetPr, mxCreateFull
integer nk, nx, nz, nh, N, Ts

! check for proper number of arguments. 
if(nrhs .ne. 17) then
	call mexErrMsgTxt('17 input variables required.')
elseif(nlhs .ne. 3) then
	call mexErrMsgTxt('3 output variables required.')
endif

! get the size of the input array
N    = mxGetM(prhs(1))
Ts   = mxGetN(prhs(3))
nk   = mxGetM(prhs(5))
nx   = mxGetM(prhs(7))
nz   = mxGetM(prhs(9))
nh   = mxGetM(prhs(10))

! create matrix for the return arguments.
plhs(1)   = mxCreateFull(N, 1, 0)
plhs(2)   = mxCreateFull(N, 1, 0)
plhs(3)   = mxCreateFull(1, Ts, 0)

call ssIEfcnP( %val(mxGetPr(plhs(1))), %val(mxGetPr(plhs(2))), %val(mxGetPr(plhs(3))),                             &
                %val(mxGetPr(prhs(1))), %val(mxGetPr(prhs(2))), %val(mxGetPr(prhs(3))), %val(mxGetPr(prhs(4))),     &
	            %val(mxGetPr(prhs(5))), nk, %val(mxGetPr(prhs(7))), nx, %val(mxGetPr(prhs(9))),                     &
		        %val(mxGetPr(prhs(10))), nh, N, Ts, %val(mxGetPr(prhs(14))), %val(mxGetPr(prhs(15))),               &
			    %val(mxGetPr(prhs(16))), %val(mxGetPr(prhs(17))), nz )

return
end

!   1    2     3               1    2    3   4    5  6   7  8   9 10  11 12  13   14     15   16    17
! [kpd, zpd, simh] = ssIEfcnP(kd0, zd0, sx, optK, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz)

subroutine ssIEfcnP(kpd, zpd, simh, kd0, zd0, sx, optK, k, nk, x, nx, z, h, nh, N, Ts, alpha, eta, rhoz, stdz, nz)

use dfwin 
use NUMERICAL_LIBRARIES
implicit none


! specify input and output variables
integer, intent(in) :: nk, nx, nz, nh, N, Ts
real*8, intent(out) :: Kpd(N), zpd(N), simh(Ts)
real*8, intent(in) :: kd0(N), zd0(N), sx(Ts), optK(nk, nh, nx, nz), k(nk), x(nx), z(nz), h(nh)
real*8, intent(in) :: alpha, eta, rhoz, stdz

! specify intermediate variables
real*8  :: kd(N), zd(N), optK_x(nk, nh, nz), optK_xh(nk, nz), optK_xhk(nz)
real*8  :: shocHalf(N/2), shocFull(N)
real*8  :: T1, lapse, pup, kptmp
integer :: t, len, up, down, j
logical*4 :: statConsole, xloc(nx), kloc(nk), hloc(nh)


! initializing
kpd = kd0
zpd = zd0

T1 = secnds(0.0)

do t = 1, Ts
    ! update current period capital and z shock
    kd       = kpd
    zd       = zpd
    
    ! log output price --- N is a scaling factor
    simh(t) = -eta*dlog(sum(dexp(sx(t) + zd)*(kd**alpha)) /N)

    ! update next period capital stock --- tough 4-D interpolation
    ! 
    ! step 1 -- linear interpolation (extrapolation) along x dimension -- result is [nk nh nz] submatrix
    xloc  = (x >= sx(t))
	len   = count(xloc)
    if (len == 0) then
        ! sx(t) is larger than all points on x grid
        optK_x = optK(:, :, nx, :) + ((sx(t) - x(nx))/(x(nx) - x(nx - 1)))*(optK(:, :, nx, :) - optK(:, :, nx - 1, :))
    else if (len == nx) then
        ! sx(t) is smaller than all points on x grid 
        optK_x = optK(:, :, 1, :) - ((x(1) - sx(t))/(x(2) - x(1)))*(optK(:, :, 2, :) - optK(:, :, 1, :))
    else
        ! sx(t) lies in the middle of x grid    
        up     = nx - len
		down   = up + 1
        pup    = (x(down) - sx(t))/(x(down) - x(up))
        optK_x = pup*optK(:, :, up, :) + (1 - pup)*optK(:, :, down, :)
    end if
    ! 
    ! step 2 -- linear interpolation along h dimension -- result is [nk nz] submatrix
    hloc    = (h >= simh(t))
	len     = count(hloc)
    if (len == nh) then
	    ! required simh is to the left to the h grid
        optK_xh = optK_x(:, 1, :)
    else if (len == 0) then
	    ! simh is larger than all the h grid
        optK_xh = optK_x(:, nh, :)
    else              
        ! distance of two adjunct grid points
        up      = nh - len
        down    = up + 1
        ! probability of right grid point gets picked is the relative distance to the left point
        pup     = (h(down) - simh(t))/(h(down) - h(up))
        ! linear interpolation along h dimension
        optK_xh = pup*optK_x(:, up, :) + (1 - pup)*optK_x(:, down, :)
    end if
    !
    ! step 3 -- point-by-point bilinear interpolation (extrapolation) along kd and zd dimension
    do j = 1, N
        ! find up and down index of kd(j) on k grid
        kloc = (k >= kd(j))
        down = nk - count(kloc) + 1
        if (down == 1) then 
            optK_xhk = optK_xh(1, :)
        else
            ! linear interpolation along k dimension (continuous space)
            up       = down - 1
            pup      = (k(down) - kd(j))/(k(down) - k(up))
            optK_xhk = pup*optK_xh(up, :) + (1 - pup)*optK_xh(down, :)
        end if
        ! linear interpolation along z dimension (continuous space)
		call interp1(kptmp, z, optK_xhk, zd(j), nz, 1, 1)
        kpd(j) = kptmp
    end do

    ! update idiosyncratic shock in continuous state space
	call drnnoa(N/2, shocHalf)
	shocFull(1 : N/2)     = shocHalf
	shocFull(N/2 + 1 : N) = -shocHalf
    zpd = max(-3.5*(stdz/sqrt(1 - rhoz**2)), min(3.5*(stdz/sqrt(1 - rhoz**2)), rhoz*zd + stdz*shocFull))
end do

return 
end


subroutine interp1(v, x, y, u, m, n, col)
!--------------------------------------------------------------------------------------------------
! Linear interpolation routine similar to interp1 with 'linear' as method parameter in Matlab
! 
! OUTPUT:
!   v - function values on non-grid points (n by col matrix)  
! 
! INPUT: 
!   x   - grid (m by one vector) 
!   y   - function defined on the grid x (m by col matrix)
!   u   - non-grid points on which y(x) is to be interpolated (n by one vector)
!   m   - length of x and y vectors
!   n   - length of u and v vectors
!   col - number of columns of v and y matrices
! 
! Four ways to pass array arguments:
! 1. Use explicit-shape arrays and pass the dimension as an argument(most efficient)
! 2. Use assumed-shape arrays and use interface to call external subroutine
! 3. Use assumed-shape arrays and make subroutine internal by using "contains"
! 4. Use assumed-shape arrays and put interface in a module then use module
!
! INTERFACE
!   SUBROUTINE binary_search(list, element, k)
!     implicit none
!     integer*4 :: k
!     real*8    :: list(:), element
!   END SUBROUTINE binary_search
! END INTERFACE
!
! This subroutine is equavilent to the following matlab call
! v = interp1(x, y, u, 'linear', 'extrap') with x (m by 1), y (m by col), u (n by 1), and v (n by col)
!------------------------------------------------------------------------------------------------------
implicit none
 
integer :: m, n, col, i, j
real*8, intent(out) :: v(n, col)
real*8, intent(in)  :: x(m), y(m, col), u(n)
real*8    :: prob

do i = 1, n
    if (u(i) < x(1))  then
	    ! extrapolation to the left
	    v(i, :) = y(1, :) - (y(2, :) - y(1, :))   * ((x(1) - u(i))/(x(2) - x(1)))
    else if (u(i) > x(m)) then
	    ! extrapolation to the right
		v(i, :) = y(m, :) + (y(m, :) - y(m-1, :)) * ((u(i) - x(m))/(x(m) - x(m-1)))
    else
	    ! interpolation
    	! find the j such that x(j) <= u(i) < x(j+1)
	    call bisection(x, u(i), m, j)
		prob    = (u(i) - x(j))/(x(j+1) - x(j))
		v(i, :) = y(j, :)*(1 - prob) + y(j+1, :)*prob
	end if 
end do 

end subroutine interp1


subroutine bisection(list, element, m, k)
implicit none
integer*4 :: m, k, first, last, half
real*8    :: list(m), element

first = 1
last  = m
do
	if (first == (last-1)) exit
	half = (first + last)/2
	if ( element < list(half) ) then
		! discard second half
		last = half
	else
		! discard first half
		first = half
	end if
end do

k = first

end subroutine bisection
