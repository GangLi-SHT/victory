module math_mod
 
  use global_parameter
  implicit none

  integer, parameter :: n_LM = 6         ! number of parameters                                                                          
  integer, parameter :: m_LM = 6         ! number of data to be fitted       
  real(dp) :: a(m_LM), y(m_LM), b_LM

  interface simpson
     module procedure Simpson_Dble, Simpson_Cmplx
  end interface

  interface
     SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,  &
                       WORK, LWORK, RWORK, INFO)
       CHARACTER(LEN = 1), intent(in) :: JOBVL, JOBVR
       INTEGER, intent(in)  :: N, LDA, LDVL, LDVR, LWORK
       INTEGER, intent(out) :: INFO
       REAL(selected_real_kind(4)), intent(out) :: RWORK(2*N)
       COMPLEX(selected_real_kind(4)), intent(inout) :: A(LDA, N)
       COMPLEX(selected_real_kind(4)), intent(out) :: W(N), VL(LDVL, N), VR(N, LDVR)
       COMPLEX(selected_real_kind(4)), intent(out) :: WORK(LWORK)
     END SUBROUTINE CGEEV

     SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,  &
                       WORK, LWORK, RWORK, INFO)
       CHARACTER(LEN = 1), intent(in) :: JOBVL, JOBVR
       INTEGER, intent(in)  :: N, LDA, LDVL, LDVR, LWORK
       INTEGER, intent(out) :: INFO
       REAL(selected_real_kind(8)), intent(out) :: RWORK(2*N)
       COMPLEX(selected_real_kind(8)), intent(inout) :: A(LDA, N)
       COMPLEX(selected_real_kind(8)), intent(out) :: W(N), VL(LDVL, N), VR(N, LDVR)
       COMPLEX(selected_real_kind(8)), intent(out) :: WORK(LWORK)
     END SUBROUTINE ZGEEV
  end interface
  
contains
  !--------------------------------------------------------------------------
  subroutine Simpson_Cmplx(Func, nL, h, solution)
    
    implicit none
    integer, intent(in) :: nL
    complex(dp), intent(in) ::  Func(nL)
    complex(dp), intent(out) :: solution
    real(dp), intent(in) :: h

    integer :: i

    if (Mod(nL, 2) == 0) then
       print*, "nL must be odd!"
       stop
    endif

    solution = (Func(1) + Func(nL))
    do i = 2, nL-1, 2
       solution = solution + 4.d0 * Func(i)
    enddo
    do i = 3, nL-2, 2
       solution = solution + 2.d0 * Func(i)
    enddo
    solution = solution * h/3.d0

  end subroutine Simpson_Cmplx

  !-----------------------------------------------------!
  subroutine Simpson_Dble(Func, nL, h, solution)

    implicit none
    integer, intent(in)   :: nL
    real(dp), intent(in)  :: Func(nL), h
    real(dp), intent(out) :: solution

    integer :: i

    if (Mod(nL, 2) == 0) then
       print*, "nL must be odd!"
       stop
    endif

    solution = (Func(1) + Func(nL))
    do i = 2, nL-1, 2
       solution = solution + 4.d0 * Func(i)
    enddo
    do i = 3, nL-2, 2
       solution = solution + 2.d0 * Func(i)
    enddo
    solution = solution * h/3.d0

  end subroutine Simpson_Dble

  !-----------------------------------------------------------------------  
  SUBROUTINE FDfit(L, Gt, dt, FD1,FD2)

    implicit none
    integer,  intent(in) :: L
    real(dp), intent(in) :: Gt(L), dt
    real(dp), intent(out) :: FD1, FD2

    FD1 = (-25.d0*gt(1) + 48.d0*gt(2) - 36.d0*gt(3) + 16.d0*gt(4) -  3.d0*gt(5) ) &
         / 12.d0 / dt
    FD2 = ( 25.d0*gt(L) - 48.d0*gt(L-1) + 36.d0*gt(L-2) - 16.d0*gt(L-3) &
         + 3.d0*gt(L-4) ) / 12.d0 / dt
    RETURN
  END SUBROUTINE FDfit

  !------------------------------------------------------------------------
  SUBROUTINE nfourier(Mtype, L, Iwmax, FD1, FD2, rindata, coutdata)
    !
    ! Purpose
    ! =======
    !   Fourier transform using the first derivative at both the first and the end point of G(tau)
    !
    ! Arguments
    ! =========
    character(len=10), intent(in) :: Mtype
    integer,  intent(in)     :: L, Iwmax
    real(dp), intent(in)     :: FD1, FD2, rindata(L+1)
    complex(dp), intent(out) :: coutdata(Iwmax)

    ! ... local vars ... 
    integer     :: i, j, k
    real(dp)    :: rincopy(L+1), a(L), b(L), c(L), d(L), u(L+1), q(L+1), xm(L+1)
    real(dp)    :: delta, three, six, p, om
    COMPLEX(dp) :: cdummy,explus,ex

    delta = Beta/L
    DO i = 1, L + 1
       rincopy(i) = rindata(i)
    ENDDO
    three = Two+One
    six = two*three

    !     spline interpolation:  the spline is given by
    !     G(tau) = a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
    !     The following formulas are taken directly from  Stoer and
    !     Bulirsch p. 102
    IF(FD1 .GT. .99E30) THEN
!       print *,'Natural spline boundary condition is appeared! Check!'
       q(1) = Zero
       u(1) = Zero
    ELSE
       q(1) = -0.5
       u(1) = (3./delta)*( ( rincopy(2)-rincopy(1) )/delta -FD1 )
    ENDIF
    DO k = 2, L
       p = q(k-1)/Two+Two
       q(k)=-One/Two/p
       u(k)=Three/delta**2*(rincopy(k+1)+rincopy(k-1)-Two*rincopy(k))
       u(k)=(u(k)-u(k-1)/Two)/p
    ENDDO

    IF(FD2 .GT. .99E30) THEN
!       print *,'Natural spline boundary condition is appeared! Check!'
       q(L+1) = Zero
       u(L+1) = Zero
    ELSE
       q(L+1) = 0.5
       u(L+1) = (3./delta)*(FD2 - (rincopy(L+1)-rincopy(L))/delta )
    ENDIF

    XM(L+1) = ( u(L+1)-q(L+1)*u(L))/(q(L+1)*q(L)+1 )

    DO k = L,1,-1
       XM(k) = q(k)*XM(k+1) + u(k)
    ENDDO

    DO j = 1, L
       a(j) = rincopy(j)
       c(j) = XM(j)/Two
       b(j) = (rincopy(j+1)-rincopy(j))/delta - (Two*XM(j)+XM(j+1))*delta/6.
       d(j) = (XM(j+1)-XM(j))/(6.*delta)
    ENDDO

    !     The Spline multiplied by the exponential can now be exlicitely
    !     integrated. The following formulas were obtained using
    !     MATHEMATICA

    DO i = 1, Iwmax
       if(Mtype == "Fermionic") then
          om = (Two*i - One)*pi/Beta
       elseif (Mtype == "Bosonic") then
          om = Two*(i-1)*Pi/Beta
       else
          print*, "Mtype is wrong in nfourier!"
          stop
       end if
       coutdata(i) = Zero
       DO j = 1, L
          cdummy = xi*om*delta*j
          explus = exp(cdummy)
          cdummy = xi*om*delta*(j-1)
          ex = exp(cdummy)
          if (abs(Om) < 1.d-6) then
             coutdata(i) = coutdata(i) + a(j)*delta + 0.5d0*b(j)*delta*delta + &
                  c(j)*delta*delta*delta/3.d0 + d(j)*delta*delta*delta*delta/4.d0
          else
             coutdata(i) = coutdata(i) + explus*(( -six* d(j) )/om**4 +  &
                  ( Two*xi*c(j) + six*delta*xi*d(j)  )/om**3 +    &
                  ( b(j)+ Two*delta*c(j)+ three*delta**2*d(j) )/om**2 +  &
                  (- xi*a(j) - delta*xi*b(j) - delta**2*xi*c(j) - &
                  delta**3*xi*d(j))/om)

             coutdata(i) = coutdata(i) + ex*(six*d(j)/om**4 - Two*xi*c(j)/om**3 &
                  -b(j)/om**2 + xi*a(j)/om)
          end if
       ENDDO
    ENDDO

  end SUBROUTINE nfourier

  !-------------------------------------------------------------------
  subroutine fftf2d(n1, n2, a, wff1, wff2)
    ! --- Arguments ---
    integer, intent(in) :: N1, N2
    complex(dp), intent(inout) :: a(N1, N2)
    real(dp), intent(in) :: wff1(4*N1+15), wff2(4*N2+15)

    ! ... local vars ...
    integer :: i
    complex(dp) :: c1(N1), c2(N2)

    do i = 1, N2
       c1(1:n1) = A(1:n1, i)
       CALL zfftf(N1, c1, WFF1)
       A(1:n1, i) = c1(1:n1)/n1
    end do
    do i = 1, N1
       c2(1:n2) = a(i, 1:n2)
       CALL zfftf(n2, c2, wff2)
       a(i,1:n2) = c2(1:n2)/n2
    end do
    return
  end subroutine fftf2d

  subroutine fftb2d(n1, n2, a, wff1, wff2)
    ! --- Arguments ---
    integer, intent(in) :: n1, n2
    complex(dp), intent(inout) :: a(N1, N2)
    real(dp), intent(in) :: wff1(4*N1+15), wff2(4*N2+15)

    ! ... local vars ...
    integer :: i
    complex(dp) :: c1(n1), c2(n2)

    do i = 1, n2
       c1(1:n1) = A(1:n1, i)
       CALL zfftb(N1, c1, WFF1)
       A(1:n1, i) = c1(1:n1)
    end do
    do i = 1, n1
       c2(1:n2) = a(i,1:n2)
       CALL zfftb(n2, c2, wff2)
       a(i,1:n2) = c2(1:n2)
    end do
    return
  end subroutine fftb2d

end module math_mod
