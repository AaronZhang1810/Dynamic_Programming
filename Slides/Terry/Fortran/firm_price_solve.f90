!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! firm_price_solve.f90
!
! Solve the example dynamic price adjustment model for
! the Structural Summer School.
!
! Stephen Terry
!
! This Version: August 2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params
implicit none

!set model parameters
double precision, parameter :: epsilon = 2; !demand elasticity
double precision, parameter :: c = 0.75 !quadratic price adjustment cost
double precision, parameter :: beta = 0.96 !discount rate
double precision, parameter :: rho = 0.85 !pers of marginal cost
double precision, parameter :: sigma = 0.05 !std dev of marginal cost innnovation

!set grid parameters
integer, parameter :: pnum = 1000 !number of price grid points
double precision, parameter :: pmin = 0.75 !price grid lb
double precision, parameter :: pmax = 3.5 !price grid ub
integer, parameter :: mnum = 35  !number of marginal cost grid points
double precision, parameter :: mnstdev = 3.0 !number of std dev's around steady state for m
integer, parameter :: statenum = mnum*pnum !total number of state grid points

!set solution parameters
double precision, parameter :: soltol = 1.0e-7 !tolerance for iterative algorithms
integer, parameter :: maxit = 1000 !max iterations for iterative algorithms

!set simulation parameters
integer, parameter :: firmnum = 5000 !# of firms in simulation
integer, parameter :: Terg = 50 !# of burn-in periods per firm
integer, parameter :: Tsim = 15 !# of periods for useful simulation
integer, parameter :: Ttot = Terg + Tsim + 1 !total # of periods for simulation
integer, parameter :: seedint = 25477 !seed for random shocks
integer, parameter :: minit = mnum/2 !initial point for m in simulation
integer, parameter :: pinit = pnum/2 !initial point for p in simulation

!some GMM setup
integer, parameter :: nparam = 1

!some globally available stuff
integer :: seeddim
double precision :: start

double precision, allocatable :: mgrid(:),pgrid(:),V(:,:),Vold(:,:),pol(:,:),polold(:,:),&
		pr_mat_m(:,:),Rmat(:,:,:),dist(:,:),distold(:,:),distp(:),distm(:),simshockm(:,:),&
		msim(:,:),psim(:,:)

integer, allocatable :: polind(:,:),polindold(:,:),seedarray(:),msimpos(:,:),psimpos(:,:)

end module params

module solve
use params
use base_lib
use omp_lib
implicit none

contains

double precision function fGMM(x)
use params
use base_lib
use omp_lib
implicit none

!declarations
integer :: ct
double precision :: x(nparam)

!set the number of threads here
!call omp_set_num_threads(4)

!$omp parallel
write(*,*) "Parallel hello!"
!$omp end parallel

start = omp_get_wtime()
write(*,*) "#############################################"
write(*,*) "Solving the dynamic price adjustment problem"
write(*,*) " "

!write some parameters
write(*,*) "##### Parameters for this solution "
write(*,*) "epsilon = ",epsilon
write(*,*) "c = ",c
write(*,*) "beta = ",beta
write(*,*) "rho = ",rho
write(*,*) "sigma = ",sigma
write(*,*) " "

!do some simple allocations
allocate(mgrid(mnum),pgrid(pnum),V(pnum,mnum),Vold(pnum,mnum),pol(pnum,mnum),&
	polold(pnum,mnum),polind(pnum,mnum),polindold(pnum,mnum),pr_mat_m(mnum,mnum),&
	Rmat(pnum,mnum,pnum),dist(pnum,mnum),distold(pnum,mnum),distp(pnum),&
	distm(mnum),simshockm(firmnum,Ttot),msimpos(firmnum,Ttot),psimpos(firmnum,Ttot),&
	msim(firmnum,Ttot),psim(firmnum,Ttot))

!discretize the states and the shock process
call discretize()

!set up the return arrays - firm payouts
call rewardsetup()

!do the VFI
call solvemodel()

!now, compute ergodic distributions
call ergdists()

!now, do the simulation
call simulation()

!write the output files
call writefiles()

fGMM = 0.0

!deallocate the necessary arrays
deallocate(mgrid,pgrid,V,Vold,pol,polold,polind,polindold,&
	pr_mat_m,Rmat,dist,distold,distp,distm,simshockm,msimpos,&
	psimpos,msim,psim)

write(*,*) "Done with function evaluation at ",omp_get_wtime()-start," secs."
write(*,*) "#############################################"

end function fGMM


subroutine simulation()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,firmct,t,msavect,pmin1ct,mct,mprimect,pct
double precision :: runval


write(*,*) "##### Doing firm simulation "

!this block seeds the random draws
call random_seed(size=seeddim)
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint+ct
end do !ct
call random_seed(put=seedarray)
!call the random draws for the productivity process
call random_number(simshockm); !U(0,1) in each cell

deallocate(seedarray)


!initialize the simulation arrays
msimpos(:,:) = 0
psimpos(:,:) = 0

!initialize the productivity and input positions according to the ergodic distribution
t = 1
msimpos(:,t) = minit
psimpos(:,t) = pinit

!persistent first
do t=2,Ttot
do firmct=1,firmnum
	mct = msimpos(firmct,t-1)
	msavect=0
	runval = 0.0
	do mprimect=1,mnum
		runval = runval+pr_mat_m(mct,mprimect)
		if ((simshockm(firmct,t)<=runval).and.(msavect==0)) then
			msavect = mprimect
		end if
	end do !mprimect

	msimpos(firmct,t) = msavect
	msim(firmct,t) = mgrid(msavect)

end do !firmct
end do !t

!now, do endogenous simulation
do firmct=1,firmnum
	do t=2,Ttot
		pmin1ct = psimpos(firmct,t-1)
		mct = msimpos(firmct,t)
		pct = polind(pmin1ct,mct)
		psimpos(firmct,t) = pct
		psim(firmct,t) = pgrid(pct)
	end do !t
end do !firmct

write(*,*) "Simulation bounds check: p lb, min(p sim), max(p sim), p ub"
write(*,*) pgrid(1),minval(psim(:,(Terg+1):(Ttot-1))),maxval(psim(:,(Terg+1):(Ttot-1))),pgrid(pnum)
write(*,*) " "

end subroutine simulation

subroutine ergdists()
use params
use base_lib
use omp_lib
implicit none

integer :: pmin1ct,pct,mct,distct,mprimect
double precision :: disterr


write(*,*) "##### Computing ergodic distribution "


!initialize distribution
distold(:,:)  = 1.0
distold = distold/sum(distold)

!now do the distributional iteration
do distct=1,maxit

	dist(:,:) = 0.0

	!loop over states
	do pmin1ct=1,pnum
		do mct=1,mnum
			if (distold(pmin1ct,mct)>0.0) then

				!extract policy
				pct = polind(pmin1ct,mct)

				!push weight forward
				do mprimect=1,mnum
					dist(pct,mprimect) = dist(pct,mprimect) + &
					 	pr_mat_m(mct,mprimect) * distold(pmin1ct,mct)
				end do !mprimect

			end if
		end do !mct
	end do !pmin1ct

	!compute error
	disterr = maxval(abs(dist-distold))

	if (mod(distct,5)==1) then
		write(*,*) "Dist iter = ",distct," dist err = ",disterr
	end if

	!if dist converged, then exit
	if (disterr<soltol) then
		exit
	end if

	!otherwise, iterate forward
	distold = dist
	distold = distold/sum(distold)

end do !distct

!now, compute marginal ergodic distributions
distp(:) = 0.0
do pct=1,pnum
	distp(pct) = sum(dist(pct,:))
end do !pct


distm(:) = 0.0
do mct=1,mnum
	distm(mct) = sum(dist(:,mct))
end do !mct

!output some diagnostics
write(*,*) "Ergodic distribution for p: wgt at plb, wgt at pub"
write(*,*) distp(1),distp(pnum)
write(*,*) " "


end subroutine ergdists

subroutine solvemodel()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,pmin1ct,mct,pct,vfct,mprimect
double precision :: RHSvec(pnum),pmin1val,mval,pval,vferr,polerr,maxshift,minshift

write(*,*) "##### Implementing VFI "

!initialize guesses for vf
V(:,:) = 0.0 !(p_{-1},m)
Vold(:,:) = 0.0
pol(:,:) = 0.0
polold(:,:) = 0.0
polindold(:,:) = 0

!start the VFI loop
do vfct=1,maxit

	!within each iteration, loop over states
	!$omp parallel private(pmin1ct,mct,pmin1val,mval,pval,RHSvec,pct,mprimect)
	!$omp do collapse(2)
	do pmin1ct=1,pnum
		do mct=1,mnum

			!extract states
			pmin1val = pgrid(pmin1ct)
			mval = mgrid(mct)
			pval = pgrid(pct)

			!initialize vector of potential choices
			RHSvec(:) = 0.0

			!consider candidate policies pct
			do pct=1,pnum

				!first, insert static payoff
				RHSvec(pct) = Rmat(pmin1ct,mct,pct)

				!now, increment with expected continuation value
				do mprimect=1,mnum
					RHSvec(pct) = RHSvec(pct) + beta * &
					 	pr_mat_m(mct,mprimect) * &
						Vold(pct,mprimect)
				end do !mprimect

			end do !pct

			!extract max value
			pct = maxloc(RHSvec,1)
		  polind(pmin1ct,mct) = pct
			pol(pmin1ct,mct) = pgrid(pct)
	   	V(pmin1ct,mct) = RHSvec(pct)

		end do !mct
	end do !pmin1ct
	!$omp end do nowait
	!$omp end parallel

	!compute errors
	vferr = maxval(abs(V-Vold))
	polerr = maxval(abs(pol-polold))

	if (mod(vfct,5)==1) then
		write(*,*) "VF iter = ",vfct," VF Err = ",vferr," Pol Err = ",polerr
	end if

	!if vf converged, then exit
	if (vferr<soltol) then
		exit
	end if

	!if not converged, then update using MQP bounds and proceed
	minshift = (beta/(1-beta))*minval(V-Vold);
  maxshift = (beta/(1-beta))*maxval(V-Vold);
	!Vold = V + (maxshift+minshift)/2;
	Vold = V
	polold = pol
	polindold = polind

end do !vfct

!do some diagnostics
write(*,*) "p lb, min(pol), max(pol), p ub"
write(*,*) pgrid(1),minval(pol),maxval(pol),pgrid(pnum)

write(*,*) " "

end subroutine solvemodel

subroutine rewardsetup()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,pmin1ct,mct,pct

write(*,*) "##### Setting up the return arrays "

!set up the manager reward function
Rmat(:,:,:) = 0.0 !(p,m,p')

do pmin1ct=1,pnum
do mct=1,mnum
do pct=1,pnum

	Rmat(pmin1ct,mct,pct) = payoff(pgrid(pct),mgrid(mct),pgrid(pmin1ct))

end do !pct
end do !mct
end do !pmin1ct

write(*,*) " "

end subroutine rewardsetup

subroutine writefiles()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,pmin1ct,mct,pct,t,firmct

open(13,file="constantvec.txt")
write(13,*) mnum,pnum,Ttot,firmnum,statenum,Terg,Tsim
close(13)

open(13,file="pr_mat_m.txt")
do ct=1,mnum
write(13,*) pr_mat_m(ct,:)
end do !ct
close(13)

open(13,file="mgrid.txt")
do ct=1,mnum
write(13,*) mgrid(ct)
end do !ct
close(13)

open(13,file="pgrid.txt")
do ct=1,pnum
write(13,*) pgrid(ct)
end do !ct
close(13)

open(13,file="V.txt")
do pmin1ct=1,pnum
	do mct=1,mnum
		write(13,*) V(pmin1ct,mct)
	end do !mct
end do !pmin1ct
close(13)

open(13,file="pol.txt")
do pmin1ct=1,pnum
	do mct=1,mnum
		write(13,*) pol(pmin1ct,mct)
	end do !mct
end do !pmin1ct
close(13)

open(13,file="dist.txt")
do pmin1ct=1,pnum
	do mct=1,mnum
		write(13,*) dist(pmin1ct,mct)
	end do !mct
end do !pmin1ct
close(13)

open(13,file="distp.txt")
do ct=1,pnum
write(13,*) distp(ct)
end do !ct
close(13)

open(13,file="distm.txt")
do ct=1,mnum
write(13,*) distm(ct)
end do !ct
close(13)

open(13,file="msimpos.txt")
do firmct=1,firmnum
	do t=1,Ttot
		write(13,*) msimpos(firmct,t)
	end do !t
end do !firmct
close(13)

open(13,file="msim.txt")
do firmct=1,firmnum
	do t=1,Ttot
		write(13,*) msim(firmct,t)
	end do !t
end do !firmct
close(13)

open(13,file="psimpos.txt")
do firmct=1,firmnum
	do t=1,Ttot
		write(13,*) psimpos(firmct,t)
	end do !t
end do !firmct
close(13)

open(13,file="psim.txt")
do firmct=1,firmnum
	do t=1,Ttot
		write(13,*) psim(firmct,t)
	end do !t
end do !firmct
close(13)

end subroutine writefiles

subroutine discretize()
use params
use base_lib
use omp_lib
implicit none

integer :: ct,mct,pct

write(*,*) "##### Discretizing the model "


!first, discretize the capital grid with log-linear spacing
call linspace(pgrid,log(pmin),log(pmax),pnum)
pgrid = exp(pgrid)

!now, discretize the persistent component eps
call idio_tauchen(mnum,rho,sigma,pr_mat_m,mgrid,mnstdev)



write(*,*) " "

end subroutine discretize

double precision function ACp(pval,pmin1val)
use params
use base_lib
use omp_lib
implicit none

double precision :: pval,pmin1val

ACp =  (c/2.0) * ( pval-pmin1val )**2.0

end function ACp

double precision function Piflow(pval,mval)
use params
use base_lib
use omp_lib
implicit none

double precision :: pval,mval

Piflow =  (pval-mval)*(pval**(-1.0*epsilon))

end function Piflow

double precision function payoff(pval,mval,pmin1val)
use params
use base_lib
use omp_lib
implicit none

double precision :: pval,mval,pmin1val,ACpval,Piflowval

ACpval = ACp(pval,pmin1val)
Piflowval = Piflow(pval,mval)
payoff = Piflowval - ACpval


end function payoff

end module solve
