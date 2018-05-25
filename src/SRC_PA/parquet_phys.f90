module Parquet_Phys
  !
  ! Purpose
  ! ========
  ! from the converged self-energy and complete vertex function, one can determine the eigenvalue of
  ! the BSE equation. When the leading eigenvalue is greater than 1, the BSE equation is clearly not converged, which
  ! signals the phase transition in the corresponding channel.  
  !
  ! Here, we use the lapack routine to diagonalize the eigen equations, which is usually enough for matrices of sizes 
  ! smaller than 4000 x 4000. Otherwise, different routine shall be adopted.  
  !
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util

contains
  !----------------------------------------------------------------------------------------------
  subroutine solve_eigen_equation
    !
    ! The eigen equation is only constucted for qw = 0. As in different channels, the leading eigenvalue can appear
    ! at different momentum, thus, it is nessary to solve the eigen equation at different momentum-q. According to our
    ! way of parallization, different momentum-q with qw=0 scatter at different nodes. The eigen equation will only be
    ! solved at these special nodes. 
    !    
 
    ! ... local vars ...
    integer       :: i, j, k, iChannel
    type(Indxmap) :: ComIdx1, ComIdx2, ComIdx3
    complex(dp)   :: dummy

    integer  :: info, LWORK
    integer, parameter :: LWMAX=10000
    real(dp) :: RWork(2*Nt)
    complex(dp) :: EigenVa(Nt), EigenVL(Nt, Nt)
    complex(dp) :: EigenVR(Nt, Nt), work(LWMAX)
    character(len=30) :: FLE

    ! ... executable ...
    do k = 1, Nb
       ComIdx1 = index_bosonic(id*Nb+k)
       if (ComIdx1%iw == 1) then   ! qw = 0
          do iChannel = 1, 3
             do j = 1, Nt
                select case (iChannel)
                case(1,2)
                   call index_operation(index_fermionic(j), Comidx1, FaddB, ComIdx2)
                   if (ComIdx2%iw > Nf .or. Comidx2%iw < 1) then
                      dummy = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1)+One) + mu - Ek(ComIdx2%ix, ComIdx2%iy))
                   else
                      dummy = Gkw(list_index(ComIdx2, Fermionic))
                   end if
                case (3,4)
                   call index_operation(Index_fermionic(j), Index_fermionic(j), MinusF, ComIdx2)    !  -k
                   call index_operation(ComIdx2, ComIdx1, FaddB, ComIdx3)                ! q-k
                   if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                      ! use the non-interacting green's function when q-k is outside the box
                      dummy = One/( xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) )
                   else
                      dummy = Gkw(list_index(ComIdx3, Fermionic))
                   end if
                end select
                do i = 1, Nt
                   select case (iChannel)
                   case (1)
                      mat(i, j, 1) =  G_d(i, j, k)*Gkw(j)*dummy/Nc/beta
                   case (2)
                      mat(i, j, 1) =  G_m(i, j, k)*Gkw(j)*dummy/Nc/beta
                   case (3)
                      mat(i, j, 1) =  Half*(G_s(i, j, k)-G_t(i, j, k))*Gkw(j)*dummy/Nc/beta
                   end select
                end do
             end do
             
             ! Query the optimal workspace.
             LWORK = -1
             call zgeev('Vectors', 'Vectors', Nt, mat(1:Nt,1:Nt,1), Nt, EigenVa, EigenVL, Nt, EigenVR, &
                  Nt, Work, LWORK, Rwork, info)
             LWORK = MIN(LWMAX, INT( WORK(1) ))
             call zgeev('Vectors', 'Vectors', Nt, mat(1:Nt,1:Nt,1), Nt, EigenVa, EigenVL, Nt, EigenVR, &
                  Nt, Work, LWORK, Rwork, info)
             
             if (info /= 0) then
                print*, "ZGEEV is wrong in calculating EigenValues!"
                stop
             end if

             select case (ichannel)
             case (1)
                write(FLE, '(a)') 'EigenVal_d.dat'
             case (2)
                write(FLE, '(a)') 'EigenVal_m.dat'
             case (3)
                write(FLE, '(a)') 'EigenVal_pp.dat'
             end select

             open(unit=132, file=FLE, status='unknown', position='append')
             do i = 1, Nt
                if (dble(EigenVa(i)) > Zero .and. dble(EigenVa(i)) < One+1.d-8) then
                   write(132, '(2i5, 5f12.6)') ComIdx1%ix, ComIdx1%iy, dble(EigenVa(i)), imag(EigenVa(i))
                end if
             end do
             close(132)             
             
          end do
       end if
    end do

  end subroutine solve_eigen_equation
  
  ! -----------------------------------------------------------------------------------------------------------
end module Parquet_Phys
