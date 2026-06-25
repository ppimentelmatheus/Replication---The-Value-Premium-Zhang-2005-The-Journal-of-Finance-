subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!-----------------------------------------------------------------------------------------------------
! Also return Zf in the list: the only difference from PANIEFCN.f90
!-----------------------------------------------------------------------------------------------------
implicit none

integer plhs(*), prhs(*)
integer nlhs, nrhs
integer mxGetM, mxGetN, mxGetPr, mxCreateFull
integer nk, nz, nx, N, Ts, nh

! check for proper number of arguments. 
if (nrhs .ne. 25) then
	call mexErrMsgTxt('25 input variables required.')
else if (nlhs .ne. 15) then
	call mexErrMsgTxt('15 output variables required.')
end if

! get the size of the input array
N    = mxGetM(prhs(1))
Ts   = mxGetN(prhs(6))
! additional shape information to be send to subroutine to facilitate variable declaration
nk   = mxGetM(prhs(5))
nh   = mxGetM(prhs(7))
nx   = mxGetM(prhs(9))
nz   = mxGetM(prhs(11))

! create matrix for the return arguments.
plhs(1)   = mxCreateFull(N, Ts, 0)
plhs(2)   = mxCreateFull(N, Ts, 0)
plhs(3)   = mxCreateFull(N, Ts, 0)
plhs(4)   = mxCreateFull(N, Ts - 1, 0)
plhs(5)   = mxCreateFull(N, Ts, 0)
plhs(6)   = mxCreateFull(1, 1, 0)
plhs(7)   = mxCreateFull(1, 1, 0)
plhs(8)   = mxCreateFull(1, 1, 0)
plhs(9)   = mxCreateFull(1, 1, 0)
plhs(10)   = mxCreateFull(1, 1, 0)
plhs(11)  = mxCreateFull(1, Ts - 1, 0)
plhs(12)  = mxCreateFull(1, Ts - 1, 0)
plhs(13)  = mxCreateFull(N, 1, 0)
plhs(14)  = mxCreateFull(N, 1, 0)
plhs(15)  = mxCreateFull(N, Ts, 0)

call panIEfcnZ( %val(mxGetPr(plhs(1))), %val(mxGetPr(plhs(2))), %val(mxGetPr(plhs(3))), %val(mxGetPr(plhs(4))),    &
                %val(mxGetPr(plhs(5))), %val(mxGetPr(plhs(6))), %val(mxGetPr(plhs(7))), %val(mxGetPr(plhs(8))),    &
                %val(mxGetPr(plhs(9))), %val(mxGetPr(plhs(10))), %val(mxGetPr(plhs(11))), %val(mxGetPr(plhs(12))), &
				%val(mxGetPr(plhs(13))), %val(mxGetPr(plhs(14))), %val(mxGetPr(plhs(15))),                         &
                %val(mxGetPr(prhs(1))), %val(mxGetPr(prhs(2))), %val(mxGetPr(prhs(3))), %val(mxGetPr(prhs(4))),    &
                %val(mxGetPr(prhs(5))), %val(mxGetPr(prhs(6))), %val(mxGetPr(prhs(7))), nh,                        &
                %val(mxGetPr(prhs(9))), nx, %val(mxGetPr(prhs(11))), N, Ts, %val(mxGetPr(prhs(14))),               & 
			    %val(mxGetPr(prhs(15))), %val(mxGetPr(prhs(16))), %val(mxGetPr(prhs(17))),                         &
			    %val(mxGetPr(prhs(18))), %val(mxGetPr(prhs(19))), %val(mxGetPr(prhs(20))),                         &  
			    %val(mxGetPr(prhs(21))), %val(mxGetPr(prhs(22))), %val(mxGetPr(prhs(23))),                         & 
			    %val(mxGetPr(prhs(24))), %val(mxGetPr(prhs(25))), nk, nz )  
 
return
end

!   1   2   3   4   5    6     7     8    9   10  11    12    13   14  15
! [Pf, Bf, Df, Rf, In,  iyr, ikr, theta, dyr, fyr, Rm, GDPg, kpd, zpd, Zf] = ...
!           1    2     3   4   5   6  7   8  9  10 11  12 13   14    15    16    17     18     19 20  21  22   23     24    25
! panIEfcn(kd0, zd0, optK, V0, k, sx, h, nh, x, nx, z, N, Ts, alpha, alp1, alp2, alp3, delta, eta, f, gP, gN, istar, rhoz, stdz)

subroutine panIEfcnZ(Pf, Bf, Df, Rf, In, iyr, ikr, theta, dyr, fyr, Rm, GDPg, kpd, zpd, Zf, &
                     kd0, zd0, optK, V0, k, sx, h, nh, x, nx, z, N, Ts,                     & 
					 alpha, alp1, alp2, alp3, delta, eta, f, gP, gN, istar, rhoz, stdz, nk, nz)

use dfwin 
use NUMERICAL_LIBRARIES
implicit none


! specify input and output variables
integer*4, intent(in) :: nk, nx, nz, N, Ts, nh
real*8, intent(out) :: Pf(N, Ts), Bf(N, Ts), Df(N, Ts), Rf(N, Ts-1), In(N, Ts), Zf(N, Ts), &
                       iyr, ikr, dyr, theta, fyr, Rm(Ts-1), GDPg(Ts-1), kpd(N), zpd(N)
real*8, intent(in) :: kd0(N), zd0(N), optK(nk, nh, nx, nz), V0(nk, nh, nx, nz), k(nk), sx(Ts), x(nx), z(nz), h(nh)
real*8, intent(in) :: alpha, alp1, alp2, alp3, gP, gN, delta, eta, f, istar, rhoz, stdz

! specify intermediate variables
real*8  :: kd(N), shocHalf(N/2), shocFull(N), Af(N, Ts), Rvf(N, Ts), tmpI(N), Vf(N, Ts)
real*8  :: V_x(nk, nh, nz), V_xh(nk, nz), V_xhk(nz), optK_x(nk, nh, nz), optK_xh(nk, nz), optK_xhk(nz)
real*8  :: T1, lapse, pup, kptmp, vtmp, simh(Ts), GDP(Ts)
integer*4 :: t, len, up, down, j
logical*4 :: statConsole, xloc(nx), kloc(nk), hloc(nh)

! initializing
kpd = kd0
zpd = zd0
! initializing the firm distribution ('d' denotes values in sampling distribution)
In  = 0.0                                         ! panel of investment
Rvf = 0.0                                         ! panel of output
Af  = 0.0                                         ! panel of adjustment cost
Df  = 0.0                                         ! panel of dividend
Vf  = 0.0                                         ! panel of firm value  
Bf  = 0.0                                         ! panel of book value (capital stock)
Zf  = 0.0                                         ! panel of idiosyncratic shock

T1 = secnds(0.0)

do t = 1, Ts
    ! update current period capital and z shock
    Bf(:, t) = kpd
    Zf(:, t) = zpd

    ! log output price --- N is a scaling factor
    simh(t) = -eta*dlog(sum(dexp(sx(t) + Zf(:, t))*(Bf(:, t)**alpha)) /N)

    ! update next period capital stock --- tough 4-D interpolation
    ! 
    ! step 1 -- linear interpolation (extrapolation) along x dimension -- result is [nk nh nz] submatrix
    xloc  = (x >= sx(t))
	len   = count(xloc)
    if (len == 0) then
        ! sx(t) is larger than all points on x grid
        optK_x = optK(:, :, nx, :) + ((sx(t) - x(nx))/(x(nx) - x(nx - 1)))*(optK(:, :, nx, :) - optK(:, :, nx - 1, :))
        V_x    = V0(:, :, nx, :) + ((sx(t) - x(nx))/(x(nx) - x(nx - 1)))*(V0(:, :, nx, :) - V0(:, :, nx - 1, :))
    else if (len == nx) then
        ! sx(t) is smaller than all points on x grid 
        optK_x = optK(:, :, 1, :) - ((x(1) - sx(t))/(x(2) - x(1)))*(optK(:, :, 2, :) - optK(:, :, 1, :))
        V_x    = V0(:, :, 1, :) - ((x(1) - sx(t))/(x(2) - x(1)))*(V0(:, :, 2, :) - V0(:, :, 1, :))
    else
        ! sx(t) lies in the middle of x grid    
        up     = nx - len
		down   = up + 1
        pup    = (x(down) - sx(t))/(x(down) - x(up))
        optK_x = pup*optK(:, :, up, :) + (1 - pup)*optK(:, :, down, :)
        V_x    = pup*V0(:, :, up, :) + (1 - pup)*V0(:, :, down, :)
    end if
    ! 
    ! step 2 -- linear interpolation along h dimension -- result is [nk nz] submatrix
    hloc    = (h >= simh(t))
	len     = count(hloc)
    if (len == nh) then
	    ! required simh is to the left to the h grid
        optK_xh = optK_x(:, 1, :)
        V_xh    = V_x(:, 1, :)
    else if (len == 0) then
	    ! simh is larger than all the h grid
        optK_xh = optK_x(:, nh, :)
        V_xh    = V_x(:, nh, :)
    else              
        ! distance of two adjunct grid points
        up      = nh - len
        down    = up + 1
        ! probability of right grid point gets picked is the relative distance to the left point
        pup     = (h(down) - simh(t))/(h(down) - h(up))
        ! linear interpolation along h dimension
        optK_xh = pup*optK_x(:, up, :) + (1 - pup)*optK_x(:, down, :)
        V_xh    = pup*V_x(:, up, :) + (1 - pup)*V_x(:, down, :)
    end if
    !
    ! step 3 -- point-by-point bilinear interpolation (extrapolation) along kd dimension
    do j = 1, N
        ! find up and down index of Bf(j, t) on k grid
        kloc = (k >= Bf(j, t))
        down = nk - count(kloc) + 1
        if (down == 1) then 
            optK_xhk = optK_xh(1, :)
            V_xhk    = V_xh(1, :)
        else
            ! linear interpolation along k dimension (continuous space)
            up       = down - 1
            pup      = (k(down) - Bf(j, t))/(k(down) - k(up))
            optK_xhk = pup*optK_xh(up, :) + (1 - pup)*optK_xh(down, :)
            V_xhk    = pup*V_xh(up, :) + (1 - pup)*V_xh(down, :)
        end if
        ! linear interpolation along z dimension (continuous space)
		call interp1(kptmp, z, optK_xhk, Zf(j, t), nz, 1, 1)
		call interp1(vtmp, z, V_xhk, Zf(j, t), nz, 1, 1)
        kpd(j)   = kptmp
		Vf(j, t) = vtmp
    end do

    ! update cross-sectional dividends
    In(:, t) = kpd - (1 - delta)*Bf(:, t)
	tmpI     = In(:, t)/Bf(:, t) - istar
    where (tmpI >= 0)
   	    Af(:, t) = (gP/2)*(tmpI**2)*Bf(:, t)
	elsewhere
	    Af(:, t) = (gN/2)*(tmpI**2)*Bf(:, t)
    end where
    ! Rvf is revenue not output
	Rvf(:, t) = dexp(sx(t) + Zf(:, t) + simh(t))*(Bf(:, t)**alpha)

    ! update idiosyncratic shock in continuous state space
	call drnnoa(N/2, shocHalf)
	shocFull(1 : N/2)     = shocHalf
	shocFull(N/2 + 1 : N) = -shocHalf
    zpd = max(-3.5*(stdz/dsqrt(1 - rhoz**2)), min(3.5*(stdz/dsqrt(1 - rhoz**2)), rhoz*Zf(:, t) + stdz*shocFull))
end do

! dividend
Df   = Rvf - In - Af - f
! important step: redefine firm value to be ex dividend to conform to empirical studies
Pf   = Vf - Df
! cross-sectional stock return via ex dividend convention
Rf   = (Pf(:, 2:Ts) + Df(:, 2:Ts)) /Pf(:, 1:Ts-1)   
! value-weighted market return
Rm   = sum(Pf(:, 1:Ts-1)*Rf, 1) /sum(Pf(:, 1:Ts-1), 1)
! GDP growth
GDP  = sum(Rvf, 1)
GDPg = GDP(2:Ts)/GDP(1:Ts-1)

! ratios
iyr   = sum(sum(In, 1)  /sum(Rvf, 1))/Ts
ikr   = sum(sum(In, 1)  /sum(Bf, 1))/Ts
theta = (sum(sum(Af, 1) /sum(Rvf, 1))/Ts) /iyr
! theta = sum(Af/In)/(Ts**2)
dyr   = sum(sum(Df, 1) /sum(Rvf, 1))/Ts
fyr   = sum(N*f/sum(Rvf, 1)) /Ts

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
