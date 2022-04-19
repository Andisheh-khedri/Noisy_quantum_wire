!-----------------------------------------------------------------------------------------------------
module commons

DOUBLE PRECISION,PARAMETER::pi=atan2(0.0d+00,-1.0d+00),w_c=100000.0d+00,D=100.0d+00
COMPLEX*16, PARAMETER::ic=cmplx(0.0d+00,1.0d+00)

end module
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!------------------------------------------------------------------
program Kane_Fischer_RG

use commons

implicit none


common nu_c,omega_c,q_c,KLL,Lnu,sfact,wprime,counter,prop

!------------------------------------------------------------------
INTEGER,PARAMETER::NEQ=1,Nsteps=8600

INTEGER:: openstatus, Flag, ISTATE, IOPT, MF,MU, ML, LRW, LIW
DOUBLE PRECISION::y,yout, RTOL, ATOL,JAC,yfun_be

DOUBLE PRECISION,allocatable,dimension(:)::RWORK
INTEGER,allocatable,dimension(:)::IWORK

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
INTEGER,parameter::ITOL=1,nrowpd=2,NEQ1=1
DOUBLE PRECISION,DIMENSION(nrowpd,NEQ1)::pd1
!------------------------------------------------------------------
!------------------------------------------------------------------
DOUBLE PRECISION::x,xout,resol,dl,KLL
INTEGER::i,kint,div_count,istar,tint
!------------------------------------------------------------------
DOUBLE PRECISION,DIMENSION(1:NEQ)::Yfun,YDfun
DOUBLE PRECISION,DIMENSION(1:Nsteps)::Yscat,YDscat,DY
!-----------------------------------------------------------------------
DOUBLE PRECISION::omega_c,nu_c,q_c
DOUBLE PRECISION::omega,nu,q,Lnu,sfact
!-----------------------------------------------------------------------
! external functions:
!-----------------------------------------------------------------------
COMPLEX*16::propagator_wq,propagator_wq_ana,propagator_w,propagator_wx
!-----------------------------------------------------------------------
REAL*4::wprime
INTEGER::counter,V0int
DOUBLE PRECISION::prop,V_0
!-----------------------------------------------------------------------------------------------------
open(unit=1,file="flow_s=1.2.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!-----------------------------------------------------------------------------------------------------
sfact=1.2d+00

KLL=1.4d+00*sin(sfact*pi/2.0d+00)*(D/w_c)**(1.0d+00-sfact)
print*,'KLL is',KLL
!------------------------------------------------------------------
Lnu=1000.0d+00/1.0d+00!
!------------------------------------------------------------------
!------------------------------------------------------------------
!----------------------------------------------------------------------
! Solving the flow equations:
!----------------------------------------------------------------------
MF=10
ML=1
MU=1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Dimension Declaration of WORK:
!----------------------------------------------------------------------
If(MF==10)then
LRW=20+16*NEQ
LIW=20
elseif(MF==21 .or. MF==22)then
LRW=22+9*NEQ+NEQ**2 
LIW=20+NEQ
elseif(MF==24 .or. MF==25)then
LRW=22+10*NEQ+(2*ML+MU)*NEQ
LIW=20+NEQ
endif

!------------------------------------------------------------------
allocate(RWORK(LRW),IWORK(LIW))
!------------------------------------------------------------------
!------------------------------------------------------------------
Lnu=100.0d+00/1.0d+00
!------------------------------------------------------------------
!------------------------------------------------------------------
Do V0int=1,10
!------------------------------------------------------------------
V_0=5.0d+00*(10.0d+00**(-2.0d+00+4.0d+00*((V0int-1.0d+00)/9.0d+00)))
!------------------------------------------------------------------
Do kint=1,1!2,2
!------------------------------------------------------------------
KLL=((2.0)*((kint-1)))*sin(sfact*pi/2.0d+00)*(D/w_c)**(1.0d+00-sfact)
!------------------------------------------------------------------
! Initialization of the solver:
!------------------------------------------------------------------
resol=1.002d+00
!------------------------------------------------------------------

Yscat(1:Nsteps)=1.0d+00*V_0

div_count=0

counter=0

Do i=1,Nsteps
!------------------------------------------------------------------
Flag=1
!------------------------------------------------------------------
IOPT=1            
RWORK(5:10)=0.0d+00 
IWORK(5:10)=0
IWORK(6)=800000 
!------------------------------------------------------------------
Yscat(1:NEQ)=1.0d+00*V_0
!------------------------------------------------------------------
RTOL=10.0d+00**(-7.0d+00) 
ATOL=10.0d+00**(-7.0d+00)
ISTATE=1
JAC=0.0d+00
!------------------------------------------------------------------

x=D/((resol)**(i-1))
xout=D/((resol)**(i))

dl=(x-xout)/x

if(i==1)then
Yscat(1)=1.0d+00*V_0
endif

CALL Fflow (NEQ, x ,Yfun, YDfun)

if(div_count==0)then
CALL DLSODE (Fflow, NEQ, Yfun, x, xout, ITOL, RTOL, ATOL, Flag,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
endif


if(i>2)then

Yscat(i+1)=Yscat(i)*(1.0d+00+dl*(1+((Yfun(1)-Yfun_be)/dl)))

else

Yscat(i+1)=Yscat(i)

endif


write(1,*)KLL,x/w_c,Yscat(i+1),(Yfun(1)-Yfun_be)/dl,dl,Yfun(1),Yfun_be,Lnu,sfact,YDfun(1),(1+((Yfun(1)-Yfun_be)/dl))&
         ,(KLL)*((x/w_c)**(-1.0d+00+sfact))*(exp(-x/w_c))*(1.0d+00/sin(pi*sfact/2.0d+00)),prop

print*,kint,tint,KLL,sfact,x,Yscat(i+1),(Yfun(1)-Yfun_be)/dl


print*,'checking:',(KLL)*((x/w_c)**(-1.0d+00+sfact))*(exp(-x/w_c))*(1.0d+00/sin(pi*sfact/2.0d+00))&
                  ,prop

Yfun_be=Yfun(1)

!------------------------------------------------------------------
end do !i
!------------------------------------------------------------------
!------------------------------------------------------------------
end do !kint
!------------------------------------------------------------------
end do !tint
!------------------------------------------------------------------
deallocate(RWORK,IWORK)
!------------------------------------------------------------------
!-----------------------------------------------------------------------
!end if ! sup
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end program Kane_Fischer_RG
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!******* Flow equations***********************
!------------------------------------------------------------------------------------------------------------
subroutine Fflow(NEQ, x, Yfun, YDfun)

use commons

implicit none

common nu_c,omega_c,q_c,KLL,Lnu,sfact,wprime,counter,prop

INTEGER::NEQ,i
DOUBLE PRECISION,DIMENSION(1:NEQ)::Yfun,YDfun
DOUBLE PRECISION::x,nu_c,omega_c,q_c,KLL,Lnu,sfact

DOUBLE PRECISION::propagator
INTEGER::counter
!---------------------------------------------------
!------------------------------------------------------------------
INTEGER*4::ier_m,neval_m,inf_m
REAL*4::bound_m,epsabs_m,epsrel_m,abserr_m,result_m,f_m,wprime,xi,xf

REAL*4::funci,funcf
DOUBLE PRECISION::prop
!-----------------------------------------------------------------

epsabs_m=0.0d+00
epsrel_m=(10.0d+00)**(-5.0d+00)
inf_m=1
bound_m=10.0d+00**(-5.0d+00)

xf=D*1000.0d+00
xi=D*(10.0d+00**(-9.0d+00))

result_m=0.0d+00

wprime=x

funci=f_m(xi)
funcf=f_m(xf)

CALL qag(f_m,xi,xf,epsabs_m,epsrel_m,6,result_m,abserr_m,neval_m,ier_m)

propagator=result_m

prop=propagator*x

counter=counter+1

YDfun(1:NEQ)=propagator

end subroutine Fflow
!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine Evaluate_propagator(x,propagator)

use commons

implicit none

common nu_c,omega_c,q_c,KLL,Lnu,sfact,wprime

!------------------------------------------------------------------
INTEGER*4::ier_m,neval_m,inf_m
REAL*4::bound_m,epsabs_m,epsrel_m,abserr_m,result_m,f_m,wprime,xi,xf

REal*4::funci,funcf
!-----------------------------------------------------------------
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION,INTENT(OUT)::propagator

DOUBLE PRECISION::nu_c,omega_c,q_c,sfact,Lnu,KLL
!-----------------------------------------------------------------
!-----------------------------------------------------------------

epsabs_m=0.0d+00
epsrel_m=(10.0d+00)**(-5.0d+00)
inf_m=1
bound_m=10.0d+00**(-5.0d+00)

xf=D*1000.0d+00
xi=D*(10.0d+00**(-9.0d+00))

result_m=0.0d+00

wprime=x

funci=f_m(xi)
funcf=f_m(xf)

CALL qag(f_m,xi,xf,epsabs_m,epsrel_m,6,result_m,abserr_m,neval_m,ier_m)

!print*,'I am rotating to Matsubara axis:'

propagator=result_m

end subroutine
!-----------------------------------------------------------------------
! Input the spectral function:
!-----------------------------------------------------------------------
!***********************************************************************
Real*4 function f_m(x)

use commons

implicit none

common nu_c,omega_c,q_c,KLL,Lnu,sfact,wprime

real*4,INTENT(in)::x

DOUBLE PRECISION::Lnu,KLL,beta,tau,sfact,propagator
REAL*4::wprime

DOUBLE PRECISION::nu_c,omega_c,q_c
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! infinite wire - Ohmic :
!----------------------------------------------------------------------
!f_m=(2.0d+00*KLL/pi)*((x/w_c)**(sfact-1.0d+00))*(1.0d+00/((x**2)+(wprime**2)))*Exp(-x/w_c)
!----------------------------------------------------------------------
! finite wire - Ohmic-class :
!----------------------------------------------------------------------
f_m=(exp(-x/w_c))*(x**(sfact-1.0d+00))*(4.0d+00/pi)*(1.0d+00/((x**2)+(wprime**2)))&
    *(KLL**2)/(1.0d+00+(KLL**2)-(1.0d+00-(KLL**2))*cos(KLL*x*Lnu))
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end function f_m
!-----------------------------------------------------------------------
