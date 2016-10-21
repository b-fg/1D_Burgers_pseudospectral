!------------------------------------------------------------------------------------------!
!                          1D Burgers equation pseudo-spectral solver
!
!      														du   1  d               d^2
!      														-- + -- -- (u*u) = nu * ---- (u)
!      														dt   2  dx              dx^2
!
!	The analytical solution considered is
!
!     u(x,t) = 4 - 2 * nu * d/dx(phi(x,t)) / phi(x,0)
!
! where
!
!     phi(x,t) = exp(-(x-4*t     ) / (4*nu*(t+1)))
!              + exp(-(x-4*t-2*pi) / (4*nu*(t+1)))
!
!	Periodic boundary conditions are considered for a [0,2*pi] domain
!
!			u(0,t) = u(2*pi,t)
!
!	The FFTW3 library is required for the FFT and IFFT transforms.
! gfortran command to compile in Ubuntu:
! 	$ gfortran burgers1d.f90 -I/usr/include -lfftw3 -o burgers1d
!
! Distributed under the GNU GENERAL PUBLIC LICENSE.
!
! Author: B. Font Garcia
! October 2016
!------------------------------------------------------------------------------------------!

module FFTW_helper
	implicit none
	include "fftw3.f"
contains
	function FFT_1d(u,N) result(uk)
		real*8, intent(in) 		:: u(0:N-1)
		integer*8, intent(in) :: N
		complex*16 						:: uk(N/2+1)
		integer*8	:: planfw

		call dfftw_plan_dft_r2c_1d(planfw,N,u,uk,FFTW_ESTIMATE)
		call dfftw_execute_dft_r2c(planfw,u,uk)
		call dfftw_destroy_plan(planfw)
		uk = uk/N ! -- scaling required
	end function FFT_1d

	function IFFT_1d(uk,N) result(u)
		complex*16, intent(in):: uk(N/2+1)
		integer*8, intent(in) :: N
		real*8								:: u(0:N-1)
		complex*16	:: uk_copy(N/2+1)
		integer*8 	:: planbw

		uk_copy = uk ! -- input is overwritten
		call dfftw_plan_dft_c2r_1d(planbw,N,uk_copy,u,FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r(planbw,uk_copy,u)
		call dfftw_destroy_plan(planbw)
	end function IFFT_1d
end module FFTW_helper

module analytical_solution
	implicit none
	real*8, parameter :: pi = 4.*ATAN(1.0d0)
contains
	subroutine solution1(u,nu,x,t,N)
		real*8, intent(in) 		:: x(0:N-1)
		real*8, intent(in)		:: nu,t
		integer*8, intent(in)	:: N
		real*8, intent(out)		:: u(0:N-1)
		real*8 		:: phi,dphi,a,b,c
		integer*8	:: j

		do j=0,N-1
				a = (x(j)-4.d0*t)
				b = (x(j)-4.d0*t-2.d0*pi)
				c = 4.d0*nu*(t+1.d0)
				phi = exp(-a**2/c)+exp(-b**2/c)
				dphi = -2.d0*a/c*exp(-a**2/c)-2.d0*b/c*exp(-b**2/c)
				u(j) = -2.d0*nu*dphi/phi+4.d0
		end do
	end subroutine solution1
end module analytical_solution

program burgers
	! External dependencies
	use FFTW_helper
	use analytical_solution
	implicit none
	! Definitions
	complex*16, parameter :: img = (0.d0,1.d0)	! sqrt(-1) value
	! Parameters
	real*8, parameter 		:: xmin = 0.					! left bound of the physical domain
	real*8, parameter 		:: xmax = 2.*PI				! right bound of the physical domain
	real*8, parameter 		:: tmax = 0.5					! max real time
	real*8, parameter 		:: dt = 0.0001				! time step size (careful!)
	integer*8, parameter 	:: N = 32							! number of grid points (should be base 2)
	real*8, parameter 		:: nu = 0.5						! viscosity
	real*8, parameter 		:: rdt = 0.5 					! at each x% of tmax. ie: rdt=50% -> 3 outputs: t= 0,0.5*tmax,tmax
	logical,parameter			:: dealiase = .false.
	! Useful variables
	integer*8		:: i,j,k,tdt,idt,it
	real*8			:: x(0:N-1),dx,alpha,XL,kx(N/2+1),t
	real*8			:: u(0:N-1),w(0:N-1),uxt(0:N-1)
	complex*16 	:: uk(N/2+1),wk(N/2+1),Ck(N/2+1),Vk(N/2+1),rhs(N/2+1)

	! Define grid
	XL = xmax-xmin 				! -- spatial length
	dx = XL/N							! -- spatial spacing
	alpha = 2.*pi/XL			! -- basic wavenumber
	do j=0,N-1
		x(j) = xmin + j*dx
	end do
	! Define wavenumbers
	do k=1,N/2+1
		kx(k) = alpha*(k-1)
	end do
	! Initialise numerical solution with analytical solution
	call solution1(uxt,nu,x,0.d0,N) ! -- compute analitical solution for t=0
	u = uxt													! -- copy analitycal solution to numerical initial solution
	! FFT u -> uk
	uk = FFT_1d(u,N)
	! Initialise time and iterations
	t = 0.; it = 0
	! Initial ouput u = u(x,0) for both analytical and numerical solutions
	open(1,status="replace",file="output.dat",action="write")
	write(1,*)'TITLE="u"'
	write(1,*)'VARIABLES="x" "u"'
	close(1)
	tdt = nint(tmax/dt)
	idt = tdt*rdt/100
	call uxt_output(uxt,x,t,N) 	! -- output analytical solution
	call uk_output(uk,x,t,N)		! -- output numerical solution

	! Update loop
	do while(t<tmax)
		! Compute convective term (pseudo-spectral)
		! -- 1. IFFT uk -> u
		! -- Dealiase
		if(dealiase) then
			do k=1,N/2+1
				if(k>=(1./3.)*N) uk(k) = 0.
			end do
		end if
		u = IFFT_1d(uk,N)
		! -- 2. Compute nonlinear term: w = u*u (Real space)
		w = u*u
		! -- 3. FFT w -> wk
		wk = FFT_1d(w,N)
		! -- 4. Compute convective term (Fourier space)
		Ck = 0.5*img*kx*wk
		! Compute the viscous term (Fourier space)
		Vk = -nu*kx*kx*uk
		! Compute RHS: RHS = viscous term - convective term
		rhs = Vk-Ck
		! Dealise (uk, wk)
		if(dealiase) then
			do k=0,N/2
				if(k>=(1./3.)*N) rhs(k) = 0.
			end do
		end if
		! Update
		uk = uk + dt*rhs
		t = t+dt
		it = it+1
		! Ouput
		if(mod(it,idt).eq.0) then
			call solution1(uxt,nu,x,t,N)
			call uxt_output(uxt,x,t,N)
			call uk_output(uk,x,t,N)
		end if
	end do

	! Quick output check
	! Dealiase

	! if(dealiase) then
	! 	do k=1,N/2+1
	! 		if(k>=(1./3.)*N) uk(k) = 0.
	! 	end do
	! end if
	! u = IFFT_1d(uk,N)
	! print*,x(0),u(0)
	! print*,x(N-1),u(N-1)
	
	write(*,*) 'Work done! No errors.'
	stop

contains
	! TECPLOT output subrutines. The 'output.dat' will contain the temporal data for the analytical and numerical solution.
	! It can be observed the numerical solution at different time steps and compared to the analytical solution at that ouput time.
	subroutine uxt_output(u,x,t,N)
		implicit none
		real*8, intent(in) 		:: u(0:N-1),x(0:N-1),t
		integer*8, intent(in) :: N
		integer :: j

		open(1, file="output.dat", status="old", position="append", action="write")
		write(1,*)
		write(1,*)'ZONE T="analytical" ,I=',N,' ,DATAPACKING=POINT, SOLUTIONTIME =',t
		do j = 0,N-1
			write(1,102) x(j),u(j)
		end do
		close(1)
		102 format(1F8.2,E18.10)
	end subroutine uxt_output

	subroutine uk_output(uk,x,t,N)
		implicit none
		complex*16, intent(in)	:: uk(N/2+1)
		real*8, intent(in)			:: x(0:N-1),t
		integer*8, intent(in)		:: N
		complex*16	:: uk_copy(N/2+1)
		integer 		:: k
		real*8			:: u(0:N-1)

		u = IFFT_1d(uk,N)

		open(1, file="output.dat", status="old", position="append", action="write")
		write(1,*)
		write(1,*)'ZONE T="numerical" ,I=',N,' ,DATAPACKING=POINT, SOLUTIONTIME =',t
		do j = 0,N-1
			write(1,102) x(j),u(j)
		end do
		close(1)
		102 format(1F8.2,E18.10)
	end subroutine uk_output

	subroutine k_output(uk,kx,N,file)
		implicit none
		complex*16, intent(in)	:: uk(N/2+1)
		real*8, intent(in) 			:: kx(N/2+1)
		integer*8, intent(in)		:: N
		character(*),intent(in)	:: file
		integer :: k

		open(1,status="replace",file=file,action="write")
		write(1,*)'TITLE="uk"'
		write(1,*)'VARIABLES="kx" "uk_real" "uk_img"'
		write(1,*)'ZONE T="',file,'" ,I=',N/2+1,' ,DATAPACKING=POINT'
		do k = 1,N/2+1
			write(1,102) kx(k),real(realpart(uk(k))),real(imagpart(uk(k)))
		end do
		close(1)
		102 format(1F8.2,2E18.10)
	end subroutine k_output
end program burgers
