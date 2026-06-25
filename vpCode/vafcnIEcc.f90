subroutine mexFunction(nlhs, plhs, nrhs, prhs)

!--------------------------------------------------------------
!                 MEX file for VFI3FCN routine
! 
! log M_{t,t+1} = log \beta + gamma (x_t - x_{t+1})
!     gamma     = gamA + gamB (x_t - xbar)
! 
! EXTENSION ON OCTOBER 3, 2002 TO COMPUTE THE VALUE OF 
!     ASSETS-IN-PLACE
!--------------------------------------------------------------

implicit none


integer plhs(*), prhs(*)
integer nlhs, nrhs
integer mxGetM, mxGetPr, mxCreateFull
integer nk, nkp, nz, nx, nh, col_hxz, col_hz, col_xz

! check for proper number of arguments. 
if(nrhs .ne. 31) then
	call mexErrMsgTxt('31 input variables required.')
elseif(nlhs .ne. 5) then
	call mexErrMsgTxt('5 output variables required.')
endif

! get the size of the input array.
nk  = mxGetM(prhs(5))
nx  = mxGetM(prhs(7))
nz  = mxGetM(prhs(11))
nh  = mxGetM(prhs(14))
nkp = mxGetM(prhs(16))
col_hxz = nx*nz*nh
col_xz  = nx*nz
col_hz  = nz*nh

! create matrix for the return arguments.
plhs(1) = mxCreateFull(nk, col_hxz, 0)
plhs(2) = mxCreateFull(nk, col_hxz, 0)
plhs(3) = mxCreateFull(nk, col_hxz, 0)
plhs(4) = mxCreateFull(nk, col_hxz, 0)
plhs(5) = mxCreateFull(nk, col_hxz, 0)

call vafcnIEcc(%val(mxGetPr(plhs(1))), %val(mxGetPr(plhs(2))), %val(mxGetPr(plhs(3))), %val(mxGetPr(plhs(4))),        &
               %val(mxGetPr(plhs(5))),                                                                                & 
               %val(mxGetPr(prhs(1))), %val(mxGetPr(prhs(2))), %val(mxGetPr(prhs(3))), %val(mxGetPr(prhs(4))),        &
	           %val(mxGetPr(prhs(5))), nk, %val(mxGetPr(prhs(7))), %val(mxGetPr(prhs(8))), nx,                        &
		       %val(mxGetPr(prhs(10))), %val(mxGetPr(prhs(11))), nz, %val(mxGetPr(prhs(13))),                         &
		       %val(mxGetPr(prhs(14))), nh, %val(mxGetPr(prhs(16))), %val(mxGetPr(prhs(17))),                         &
		       %val(mxGetPr(prhs(18))), %val(mxGetPr(prhs(19))), %val(mxGetPr(prhs(20))), %val(mxGetPr(prhs(21))),    &
		       %val(mxGetPr(prhs(22))), %val(mxGetPr(prhs(23))), %val(mxGetPr(prhs(24))), %val(mxGetPr(prhs(25))),    &
		       %val(mxGetPr(prhs(26))), %val(mxGetPr(prhs(27))), %val(mxGetPr(prhs(28))), %val(mxGetPr(prhs(29))),    &
               %val(mxGetPr(prhs(30))), %val(mxGetPr(prhs(31))), nkp, col_hxz, col_xz, col_hz)

return
end

!-----------------------------------------------------------------------------------------------------------------
!    1   2  3   4   5                1     2     3    4   5  6   7   8    9   10 11  12  13 14  15  16  
! [optK, V, Va, I, div] = vafcnIEcc(alp1, alp2, alp3, V0, k, nk, x, xbar, nx, Qx, z, nz, Qz, h, nh, kp, ...
!                       17    18     19  20   21    22   23  24   25     26    27     28     29     30      31 
!                     alpha, beta, delta, f, gamA, gamB, gP, gN, istar, kmin, kmtrx, ksubm, hmtrx, xmtrx, zmtrx)
!-----------------------------------------------------------------------------------------------------------------

subroutine vafcnIEcc(optK, V, Va, I, div,                                                                       &
                     alp1, alp2, alp3, V0, k, nk, x, xbar, nx, Qx, z, nz, Qz, h, nh, kp,                        &
		             alpha, beta, delta, f, gamA, gamB, gP, gN, istar, kmin, kmtrx, ksubm, hmtrx, xmtrx, zmtrx, &
				     nkp, col_hxz, col_xz, col_hz)

use dfwin 
implicit none 

! specify input and output variables
integer, intent(in) :: nk, nkp, nx, nz, nh, col_hxz, col_xz, col_hz
real*8, intent(out) :: V(nk, col_hxz), Va(nk, col_hxz), optK(nk, col_hxz), I(nk, col_hxz), div(nk, col_hxz)
real*8, intent(in) :: V0(nk, col_hxz), k(nk), kp(nkp), x(nx), z(nz), Qx(nx, nx), Qz(nz, nz), h(nh)
real*8, intent(in) :: alp1, alp2, alp3, xbar, kmin, alpha, gP, gN, beta, delta, gamA, gamB, f, istar
real*8, intent(in) :: kmtrx(nk, col_hxz), ksubm(nk, col_hz), zmtrx(nk, col_hxz), xmtrx(nk, col_hxz), hmtrx(nk, col_hxz)

! specify intermediate variables
real*8  :: Res(nk, col_hxz), Obj(nk, col_hxz), ObjA(nk, col_hxz), optKold(nk, col_hxz), Vold(nk, col_hxz), & 
           VAold(nk, col_hxz), tmpEMV(nkp, col_hz), tmpI(nkp), &
           tmpObj(nkp, col_hz), tmpA(nk, col_hxz), tmpQ(nx*nh, nh), detM(nx), stoM(nx), g(nkp), tmpInd(nh, nz)
real*8  :: Qh(nh, nh, nx), Qxh(nx*nh, nx*nh), Qzxh(col_hxz, col_hxz)
real*8  :: hp, d(nh), errK, errV, errVa, T1, lapse
integer :: ix, ih, iter, optJ(col_hz), ik, iz, ind(nh, col_xz), subindex(nx, col_hz)
logical*4 :: statConsole

! construct the transition matrix for kh --- there are nx number of these transition matrix: 3-d
Qh    = 0.0
do ix = 1, nx
    do ih = 1, nh
        ! compute the predicted next period kh
        hp = alp1 + alp2*h(ih) + alp3*(x(ix) - xbar)
        ! construct transition probability vector
        d  = abs(h - hp) + 1D-32
        Qh(:, ih, ix) = (1/d)/sum(1/d)
    end do
end do

! construct the compound transition matrix over (z x h) space
! compound the (x h) space
Qxh   = 0.0
do ix = 1, nx
    call kron(tmpQ, Qx(:, ix), Qh(:, :, ix), nx, 1, nh, nh)
    Qxh(:, (ix - 1)*nh + 1 : ix*nh) = tmpQ
end do
! compound the (z x h) space: h changes the faster, followed by x, and z changes the slowest
call kron(Qzxh, Qz, Qxh, nz, nz, nx*nh, nx*nh)

! available funds for the firm
Res = dexp(xmtrx + zmtrx + hmtrx)*(kmtrx**alpha) + (1 - delta)*kmtrx - f

! initializing 
Obj     = 0.0
ObjA    = 0.0
optK    = 0.0
optKold = optK + 1.0
Vold    = V0
VAold   = V0
! Some Intermediate Variables Used in Stochastic Discount Factor
detM    = beta*dexp((gamA - gamB*xbar)*x + gamB*x**2)
stoM    = -(gamA - gamB*xbar + gamB*x)

! Intermediate index vector to facilitate submatrix extracting 
ind = reshape((/1 : col_hxz : 1/), (/nh, col_xz/))
do ix = 1, nx
    tmpInd = ind(:, ix : col_xz : nx)
	do iz = 1, nz
	    subindex(ix, (iz - 1)*nh + 1 : iz*nh) = tmpInd(:, iz)
    end do
end do

! start iterations
errK  = 1.0
errV  = 1.0
iter  = 0

T1 = secnds(0.0)

do
if (errV <= 1D-3 .AND. errK <= 1D-8) then
	exit
else
	iter = iter + 1
    do ix = 1, nx
        ! next period value function by linear interpolation: nkp by nz*nh matrix
        call interp1(tmpEMV, k, detM(ix)*(matmul(dexp(stoM(ix)*xmtrx)*Vold, Qzxh(:, subindex(ix, :)))) - ksubm, kp, &
		             nk, nkp, col_hz)
        ! maximize the right-hand size of Bellman equation on EACH grid point of capital stock
        do ik = 1, nk
            ! with istar tmpI is no longer investment but a linear transformation of that
            tmpI   = (kp - (1.0 - delta)*k(ik))/k(ik) - istar
    		where (tmpI >= 0.0)
	    	    g  = gP
		    elsewhere
		        g  = gN
		    end where
		    tmpObj = tmpEMV - spread((g/2.0)*(tmpI**2)*k(ik), 2, col_hz)
            ! direct discrete maximization
            Obj(ik, subindex(ix, :))  = maxval(tmpObj, 1)
		    optJ                      = maxloc(tmpObj, 1)
            optK(ik, subindex(ix, :)) = kp(optJ)
            ! next period value of assets-in-place by linear interpolation: 1 by nz*nh matrix
            call interp1(ObjA(ik, subindex(ix, :)), k, detM(ix)*(matmul(dexp(stoM(ix)*xmtrx)*VAold, & 
			             Qzxh(:, subindex(ix, :)))) - ksubm, (1.0 - delta)*k(ik), nk, 1, col_hz)
        end do
    end do
    ! update value function and impose limited liability condition
    V  = max(Res + Obj, kmin/10)
    Va = min(V, max(Res + ObjA, kmin/10))

    ! convergence criterion
    errK  = maxval(abs(optK - optKold))
    errV  = maxval(abs(V - Vold))
    errVa = maxval(abs(Va - VAold))
    ! revise Initial Guess
    Vold    = V
    VAold   = Va
    optKold = optK

    ! visual
	if (modulo(iter, 10) == 0) then         
		lapse = secnds(T1)      	
		statConsole = AllocConsole()
		print "(a, f10.7, a, f10.7, a, f8.1, a)", " errVa:", errVa, "   errK:", errK, "   Time:", lapse, "s"
	end if
end if
end do

! visual check on errors
lapse = secnds(T1)      	
statConsole = AllocConsole()
print "(a, f10.7, a, f10.7, a, f8.1, a)", " errVa:", errVa, "   errK:", errK, "   Time:", lapse, "s"

! optimal investment and dividend  
I    = optK - (1.0 - delta)*kmtrx
tmpA = I/kmtrx - istar
where (tmpA >= 0)
    div = Res - optK - (gP/2.0)*(tmpA**2)*kmtrx
elsewhere
    div = Res - optK - (gN/2.0)*(tmpA**2)*kmtrx  
end where

return 
end


subroutine interp1(v, x, y, u, m, n, col)
!-------------------------------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------
! find index k in list such that (list(k) <= element < list(k+1)
!--------------------------------------------------------------------------------
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


subroutine kron(K, A, B, rowA, colA, rowB, colB)
!--------------------------------------------------------------------------------
! Perform K = kron(A, B); translated directly from kron.m 
! 
! OUTPUT:
!   K -- rowA*rowB by colA*colB matrix
! 
! INPUT:
!   A -- rowA by colA matrix
!   B -- rowB by colB matrix
!   rowA, colA, rowB, colB -- integers containing shape information
!--------------------------------------------------------------------------------
implicit none

integer, intent(in) :: rowA, colA, rowB, colB
real*8, intent(in) :: A(rowA, colA), B(rowB, colB)
real*8, intent(out) :: K(rowA*rowB, colA*colB)

integer :: t1(rowA*rowB), t2(colA*colB), i, ia(rowA*rowB), ja(colA*colB), ib(rowA*rowB), jb(colA*colB)

t1 = (/ (i, i = 0, (rowA*rowB - 1)) /)
ia = int(t1/rowB) + 1
ib = mod(t1, rowB) + 1
t2 = (/ (i, i = 0, (colA*colB - 1)) /)
ja = int(t2/colB) + 1
jb = mod(t2, colB) + 1
K  = A(ia, ja)*B(ib, jb)

end subroutine kron
