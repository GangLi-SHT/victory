module Parquet_Phys
  !
  ! Purpose
  ! ========
  ! from the converged self-energy and complete vertex function, one can determine the eigenvalue of
  ! the BSE equation. When the leading eigenvalue is greater than 1, the BSE equation is clearly not converged, which
  ! signals the phase transition in the corresponding channel.  
  !
  ! Here, we use the lapack routine to diagonalize the eigen equations, which is usually enough for matrix of size 
  ! smaller than 4000 x 4000. Otherwise, different routine shall be adopted.  
  !
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util
  use parquet_kernel
  
  contains
  !----------------------------------------------------------------------------------------------
  subroutine solve_eigen_equation(ite, converged)
    !
    ! The eigen equation is only constucted for qw = 0. As in different channels, the leading eigenvalue can appear
    ! at different momentum, thus, it is nessary to solve the eigen equation at different momentum-q. According to our
    ! way of parallization, different momentum-q with qw=0 scatter at different nodes. The eigen equation will only be
    ! solved at these special nodes. 
    !    
    integer, intent(in)    :: ite
    logical, intent(in)    :: Converged
    ! ... local vars ...
    integer       :: i, j, k, iChannel, ix2, iy2, iw2, idx2, iw, idx, ix, iy, n_aux, ix_k, iy_k
    integer       :: jx, jy, jw, kqx, kqy, idx_q
    type(Indxmap) :: map_q, map_jq, map_qmk, map_j, map_mj, map_mi, map_aux, map_i, map_iq, map_kkq
    complex(dp)   :: dummy, aux, susc_kernel, aux_m, susc_kernel_m
    complex(dp)   :: dummy1D(Nt), dummy1D_pp(Nt)
    ! real(dp) :: kx, ky, dkx, dky
    integer  :: info, LWORK
    integer, parameter :: LWMAX=50000  ! this value has to be big enough
    real(dp) :: RWork(2*Nt) !, gsym_p(Nx), gsym_dxy(Nx,Ny),  gsym_dx2y2(Nx,Ny) !obsolete in version 1.1
    ! real(dp) :: norm_dxy, norm_dx2y2, norm_p
    !complex(dp) :: Ps, Pdxy, Pdx2y2, Pp    !obsolete in version 1.1
    !complex(dp) :: Susc_s, Susc_dxy, Susc_dx2y2, Susc_p  
    complex(dp) :: susc_Q_m, susc_Q_pp, susc_Q_d0, susc_Q_m0, susc_Q_pp0, GG, GG_pp
    complex(dp) :: susc_Q_d, susc_Q_d_kernel, susc_Q_d_FFT
    complex(dp) :: susc_Q_m_kernel, susc_Q_m_FFT   
    complex(dp) :: EigenVa(Nt), EigenVL(Nt, Nt)
    complex(dp) :: EigenVR(Nt, Nt), work(LWMAX)
    character(len=30) :: FLE, FLE1, FLE2, FLE3, str1, str2, str3
    
    complex(dp), allocatable :: dummy3D_1(:, :)
    complex(dp), allocatable :: GG_out(:, :, :)
  
  if (Eigen.or.Converged.or.(ite==ite_max-1)) then
    
    ! ... executable ...
    do k = 1, Nb
      idx_q = id*Nb + k
      map_q = index_bosonic_IBZ(idx_q)
      
      if (map_q%iw == 1) then ! qw = 0
        !do iChannel = 1, 3
                      
         ! select case (iChannel)
          !case(1,2)
              
            GG = Zero_c
            do j = 1, Nt
            
                map_j = index_fermionic(j)
            
                call index_FaddB(map_j, map_q, map_jq) ! k+q
                ix = map_jq%ix
                iy = map_jq%iy
                
                dummy = Zero
                
!                 if (map_jq%iw > Nf) then 
!                   do ix2 = 1,Ngrain
!                     do iy2 = 1,Ngrain
!                       dummy = dummy + One/(xi*Pi/beta*(Two*(map_jq%iw-Nf/2-1)+One) + mu - Ek_grain(ix,ix2,iy, iy2) - &
!                       Sigma_H(ix, iy))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))
!                     end do
!                   end do
!                 else if (map_jq%iw < 1) then
!                   do ix2 = 1,Ngrain
!                     do iy2 = 1,Ngrain
!                       dummy = dummy + One/(xi*Pi/beta*(Two*(map_jq%iw-Nf/2-1)+One) + mu - Ek_grain(ix,ix2,iy, iy2) - &
!                       conjg(Sigma_H(ix, iy)))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))
!                     end do
!                   end do                
!                 else
                  do ix2 = 1,Ngrain
                    do iy2 = 1,Ngrain                   
                      dummy = dummy + One/(xi*Pi/beta*(Two*(map_jq%iw-Nf/2-1)+One) + mu - Ek_grain(ix,ix2,iy, iy2)- &
                      Sigma(list_index_F(map_jq)))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))!Gkw(list_index_F(map_jq))*Gkw(j)
                    end do
                  end do
!                end if
                
                dummy1D(j)=dummy/Ngrain/Ngrain
                GG = GG + dummy/Ngrain/Ngrain
                
            end do  
                
          
          !case (3)

            GG_pp = Zero_c
            do j = 1, Nt
            
                map_j = index_fermionic(j)
                
                dummy = Zero
                call index_minusF(map_j,  map_mj)    !  -k
                call index_FaddB(map_mj, map_q, map_qmk)                ! q-k
                ix = map_qmk%ix
                iy = map_qmk%iy
                
                dummy = Zero
                
!                 if (map_qmk%iw > Nf) then
!                   do ix2 = 1,Ngrain
!                     do iy2 = 1,Ngrain 
!                       ! use the non-interacting green's function when q-k is outside the box
!                       dummy = dummy + One/( xi*Pi/beta*(Two*(map_qmk%iw-Nf/2-1) + One) + mu - Ek_grain(ix,Ngrain-ix2+1, iy,Ngrain-iy2+1) -&
!                       Sigma_H(ix, iy))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))
!                     end do
!                   end do
!                 else if (map_qmk%iw < 1) then
!                   do ix2 = 1,Ngrain
!                     do iy2 = 1,Ngrain 
!                       ! use the non-interacting green's function when q-k is outside the box
!                       dummy = dummy + One/( xi*Pi/beta*(Two*(map_qmk%iw-Nf/2-1) + One) + mu - Ek_grain(ix,Ngrain-ix2+1, iy,Ngrain-iy2+1) -&
!                       conjg(Sigma_H(ix, iy)))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))
!                     end do
!                   end do
!                 else
                  do ix2 = 1,Ngrain
                    do iy2 = 1,Ngrain 
                      dummy = dummy + One/( xi*Pi/beta*(Two*(map_qmk%iw-Nf/2-1) + One) + mu - Ek_grain(ix,Ngrain-ix2+1, iy,Ngrain-iy2+1) - &
                      Sigma(list_index_F(map_qmk)))/( xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek_grain(map_j%ix,ix2, map_j%iy,iy2) -Sigma(j))
                    end do
                  end do
!                end if
                
                dummy1D_pp(j)=dummy/Ngrain/Ngrain
                GG_pp = GG_pp + dummy/Ngrain/Ngrain
            
            end do     
          !end select
          
          
          
          ix_k=map_q%ix
          iy_k=map_q%iy
          
          write(str1, '(I0.3)') ite
          write(str2, '(I0.1)') ix_k
          write(str3, '(I0.1)') iy_k
          
          !select case (ichannel)
           ! case (1)
              write(FLE1, '(a)') 'Susc_Q_d_ite_'//trim(str2)//'_'//trim(str3)//'.dat'
           ! case (2)
              write(FLE2, '(a)') 'Susc_Q_m_ite_'//trim(str2)//'_'//trim(str3)//'.dat'
           ! case (3)
              write(FLE3, '(a)') 'Susc_Q_pp_ite_'//trim(str2)//'_'//trim(str3)//'.dat'
          !end select
          
          
          !select case (iChannel)
           ! case(1)
              !susc_Q_d0 = Zero_c
              !do i = 1,Nt
                !susc_Q_d0 = susc_Q_d0-dummy1D(i)
              !end do
              susc_Q_d0 = -GG/beta/Nc
              susc_Q_m0 = -GG/beta/Nc
              
              susc_Q_d = Zero_c
              susc_Q_m = Zero_c
              
              do i = 1,Nt
                do j = 1,Nt                  
                  susc_Q_d = susc_Q_d-dummy1D(i)*(F_d(i,j,k))*dummy1D(j)
                  susc_Q_m = susc_Q_m-dummy1D(i)*(F_m(i,j,k))*dummy1D(j)
                end do
              end do
      
              
! kernel correction to susceptibility
              susc_kernel = Zero_c
              susc_kernel_m = Zero_c

              if (.not.allocated(GG_out)) allocate(GG_out(Nx,Ny,(-f_range*Nf+1):((f_range+1)*Nf)))
               
              do jw = -f_range*Nf+1, (f_range+1)*Nf !k'
                do jx = 1, Nx
                  do jy = 1, Ny 
                    map_j = indxmap(jx, jy, jw) !k'
                    if (jw > Nf) then
                      
                      call index_FaddB(map_j, map_q, map_jq) ! k'+q
                      kqx = map_jq%ix
                      kqy = map_jq%iy
                
                      dummy = Zero_c
                
                      !if (map_jq%iw > Nf .or. map_jq%iw < 1) then (always true here)
                      do ix2 = 1,Ngrain
                        do iy2 = 1,Ngrain
                          dummy = dummy + One/(xi*Pi/beta*(Two*(map_jq%iw-Nf/2-1)+One) + mu - Ek_grain(kqx,ix2,kqy, iy2) - &
                          Sigma_H(kqx, kqy))/( xi*Pi/beta*(Two*(jw-Nf/2-1) + One) + mu - Ek_grain(jx,ix2, jy,iy2) -Sigma_H(jx,jy))
                        end do
                      end do                      
                      GG_out(jx,jy,jw) = dummy/Ngrain/Ngrain
                      
                    else if (jw < 1) then
                      
                      call index_FaddB(map_j, map_q, map_jq) ! k'+q
                      kqx = map_jq%ix
                      kqy = map_jq%iy
                
                      dummy = Zero_c
                
                      !if (map_jq%iw > Nf .or. map_jq%iw < 1) then (always true here)
                      do ix2 = 1,Ngrain
                        do iy2 = 1,Ngrain
                          dummy = dummy + One/(xi*Pi/beta*(Two*(map_jq%iw-Nf/2-1)+One) + mu - Ek_grain(kqx,ix2,kqy, iy2) - &
                          conjg(Sigma_H(kqx, kqy)))/( xi*Pi/beta*(Two*(jw-Nf/2-1) + One) + mu - Ek_grain(jx,ix2, jy,iy2) -conjg(Sigma_H(jx,jy)))
                        end do
                      end do                      
                      GG_out(jx,jy,jw) = dummy/Ngrain/Ngrain
                      
                    else  
                      j = list_index_F(map_j)
                      GG_out(jx,jy,jw) = dummy1D(j)
                    end if
                    
                  end do
                end do
              end do
              
              do jw = -f_range*Nf+1, (f_range+1)*Nf !k'
                if (jw > Nf .or. jw < 1) then
                  !do jw = -f_range*Nf+1, (f_range+1)*Nf !k'
                  do jx = 1, Nx
                    do jy = 1, Ny 
                      map_j = indxmap(jx, jy, jw) !k'
                      call index_minusF(map_j,map_mj)  ! -k'
                      call index_FaddB(map_j, map_q, map_jq) ! k'+q
                      
                      do iw = -f_range*Nf+1, (f_range+1)*Nf !k
                        do ix = 1, Nx
                          do iy = 1, Ny              
                            map_i = indxmap(ix, iy, iw) !k
                            !i = list_index_F(map_i)
                            
                            ! now we need to evaluate kernels at different index combinations
                            aux = Zero_c
                            aux_m = Zero_c
                            
                            ! (1) the most direct contribution 
                            aux = aux + kernel('d',map_i, map_j, map_q) 
                            aux_m = aux_m + kernel('m',map_i, map_j, map_q)
                            
                            ! (2) from magnetic and density channels
                            call index_FaddB(map_i, map_q, map_iq) ! k+q
                            call index_minusF(map_i,map_mi)  ! -k
                            call index_FaddF(map_j,map_mi,map_aux) !k'-k
                            
                            aux = aux -0.5d0*kernel('d',map_i, map_iq, map_aux) - 1.5d0*kernel('m',map_i, map_iq, map_aux)
                            
                            aux_m = aux_m -0.5d0*kernel('d',map_i, map_iq, map_aux) + 0.5d0*kernel('m',map_i, map_iq, map_aux)                            
                            
                            ! (3) from singlet and triplet channels                            
                            call index_FaddB(map_iq, map_j, map_kkq) ! k+k'+q
                            
                            aux = aux +0.5d0*kernel('s',map_i, map_j, map_kkq) + 1.5d0*kernel('t',map_i, map_j, map_kkq)
                            
                            aux_m = aux_m -0.5d0*kernel('s',map_i, map_j, map_kkq) + 0.5d0*kernel('t',map_i, map_j, map_kkq)
                            
                            ! now with k exchanged with k' for iw inside the box
                            if (iw>=1.and.iw <=Nf) then
                            
                              aux = aux  + kernel('d',map_j, map_i, map_q) !(1)
                              
                              aux_m = aux_m  + kernel('m',map_j, map_i, map_q)
                            
                              call index_FaddF(map_i,map_mj,map_aux) !k-k'
                              aux = aux -0.5d0*kernel('d',map_j, map_jq, map_aux) - 1.5d0*kernel('m',map_j, map_jq, map_aux) !(2)
                              
                              aux_m = aux_m -0.5d0*kernel('d',map_j, map_jq, map_aux) + 0.5d0*kernel('m',map_j, map_jq, map_aux) !(2)

                              aux = aux +0.5d0*kernel('s',map_j, map_i, map_aux) + 1.5d0*kernel('t',map_j, map_i, map_kkq) !(3)
                              
                              aux_m = aux_m -0.5d0*kernel('s',map_j, map_i, map_aux) + 0.5d0*kernel('t',map_j, map_i, map_kkq) !(3)
                            end if  
                            
                            susc_kernel = susc_kernel - GG_out(ix,iy,iw)*aux*GG_out(jx,jy,jw)
                            susc_kernel_m = susc_kernel_m - GG_out(ix,iy,iw)*aux_m*GG_out(jx,jy,jw)
                            
                          end do !iy         
                        end do !ix
                      end do !iw  
                    end do !jy
                  end do !jx
                end if !jw
              end do  !jw
              
              if (allocated(GG_out)) deallocate(GG_out)
              
! FFT correction to susceptibility (this one might not work that well...)

              susc_Q_d_FFT = -xU*Chi0_ph(idx_q)*Chi0_ph(idx_q)+xU*GG*GG/beta/beta/Nc/Nc
              susc_Q_m_FFT = -susc_Q_d_FFT
                
              susc_Q_d =  susc_Q_d/beta/beta/Nc/Nc
              susc_Q_m =  susc_Q_m/beta/beta/Nc/Nc
              
              susc_Q_d_kernel =  susc_Q_d + susc_kernel/beta/beta/Nc/Nc
              susc_Q_d_FFT = susc_Q_d_FFT + susc_Q_d_kernel
              
              susc_Q_m_kernel =  susc_Q_m + susc_kernel_m/beta/beta/Nc/Nc
              susc_Q_m_FFT = susc_Q_m_FFT + susc_Q_m_kernel
              
              open(unit=23, file=FLE1, status='unknown', position='append')
              
              write(23, '( 3i5, 8f15.6)') ite, ix_k, iy_k, dble(susc_Q_d), aimag(susc_Q_d), dble(susc_Q_d_kernel), aimag(susc_Q_d_kernel), &
               dble(susc_Q_d_FFT), aimag(susc_Q_d_FFT), dble(susc_Q_d0), aimag(susc_Q_d0)
              
              close(23) 
              
            !case(2)
              !susc_Q_m0 = Zero_c
              !do i = 1,Nt
                !susc_Q_m0 = susc_Q_m0-dummy1D(i)
              !end do
             ! susc_Q_m0 = -GG/beta/Nc
              
              !susc_Q_m = Zero_c               
              !do i = 1,Nt
                !do j = 1,Nt
                  !susc_Q_m=susc_Q_m-dummy1D(i)*(F_m(i,j,k))*dummy1D(j)
                !end do
              !end do
              
              !susc_Q_m = susc_Q_m/beta/beta/Nc/Nc    
              
              open(unit=24, file=FLE2, status='unknown', position='append')
              
              write(24, '( 3i5, 8f15.6)') ite, ix_k, iy_k, dble(susc_Q_m), aimag(susc_Q_m), dble(susc_Q_m_kernel), aimag(susc_Q_m_kernel), &
               dble(susc_Q_m_FFT), aimag(susc_Q_m_FFT), dble(susc_Q_m0), aimag(susc_Q_m0)      
              
              close(24)     
              
            !case(3)
              !susc_Q_pp0 = Zero_c
              !do i = 1,Nt
                !susc_Q_pp0 = susc_Q_pp0-dummy1D(i)
              !end do
              susc_Q_pp0 = -GG_pp/beta/Nc
              
              susc_Q_pp = Zero_c
              
              do i = 1,Nt
                do j = 1,Nt
                  susc_Q_pp=susc_Q_pp-0.5d0*dummy1D_pp(i)*(F_t(i,j,k)-F_s(i,j,k))*dummy1D_pp(j)
                end do
              end do
              
              susc_Q_pp = susc_Q_pp/beta/beta/Nc/Nc    
              
              open(unit=25, file=FLE3, status='unknown', position='append')
              
              write(25, '( 3i5, 4f15.6)') ite, ix_k, iy_k, dble(susc_Q_pp), aimag(susc_Q_pp), dble(susc_Q_pp0), aimag(susc_Q_pp0)      
              
              close(25)  
              
          !end select
          do iChannel = 1, 3
          
          if ((((map_q%ix == 1).AND.(iChannel == 3)).OR.((map_q%iy == Nx_IBZ).AND.(iChannel == 2))).OR.Converged) then         
            
            do j = 1, Nt
              do i = 1, Nt
                select case (iChannel)
                  case (1)
                    mat(i, j, 1) =  G_d(i, j, k)*dummy1D(j)/Nc/beta
                  case (2)
                    mat(i, j, 1) =  G_m(i, j, k)*dummy1D(j)/Nc/beta
                  case (3)
                    mat(i, j, 1) =  0.5d0*(G_t(i, j, k)-G_s(i, j, k))*dummy1D_pp(j)/Nc/beta
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
                write(FLE2, '(a)') 'EigenVal_ite_d.dat'
              case (2)
                write(FLE2, '(a)') 'EigenVal_ite_m.dat'
              case (3)
                write(FLE2, '(a)') 'EigenVal_ite_minus.dat'
            end select
            
            
            if ((map_q%ix == 1).AND.(ichannel==3)) then
              if (.NOT. allocated(dummy3D_1)) allocate(dummy3D_1(Nt, Nt))
              
              do i = 1,Nt
                do j = 1,Nt
                  dummy3D_1(i,j)=-0.5d0*dummy1D_pp(i)*(F_t(i,j,k)-F_s(i,j,k))*dummy1D_pp(j)
                end do
              end do
              
              open(unit=123, file='Susc_pp.dat', status='unknown')
              
              do ix2 = 1,Nx
                do iy2= 1,Ny
                  do iw2 = 1,Nf
                    idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
                    do ix = 1,Nx
                      do iy= 1,Ny
                        do iw = 1,Nf
                          idx = ((ix-1)*Ny+iy-1)*Nf + iw       
                          
                          write(123, '( 6i5, 2f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(dummy3D_1(idx2,idx)), aimag(dummy3D_1(idx2,idx))
                        end do
                      end do
                    end do
                  end do
                end do
              end do             
              
              close(123)  
              
              call output_2k(F_d, F_m, k, 'F_ph_0.dat', 127)
!               open(unit=127, file='F_ph_0.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(127, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(F_d(idx2,idx,k)), aimag(F_d(idx2,idx,k)), dble(F_m(idx2,idx,k)), aimag(F_m(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(127) 
              
              call output_2k(F_s, F_t, k, 'F_pp_0.dat', 128)
!               open(unit=128, file='F_pp_0.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(128, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(F_s(idx2,idx,k)), aimag(F_s(idx2,idx,k)), dble(F_t(idx2,idx,k)), aimag(F_t(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(128) 
              
              call output_2k(G_d, G_m, k, 'Gamma_ph_0.dat', 129)
!               open(unit=129, file='Gamma_ph_0.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(129, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(G_d(idx2,idx,k)), aimag(G_d(idx2,idx,k)), dble(G_m(idx2,idx,k)), aimag(G_m(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do
!               
!               close(129) 
              
              call output_2k(G_s, G_t, k, 'Gamma_pp_0.dat', 122)
!               open(unit=122, file='Gamma_pp_0.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(122, '(6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1,&
!                             dble(0.5d0*(G_t(idx2, idx, k)-G_s(idx2, idx, k))), aimag(0.5d0*(G_t(idx2, idx, k)-G_s(idx2, idx, k))), dble(G_s(idx2, idx, k)), aimag(G_s(idx2, idx, k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(122)  
              
              open(unit=133, file='Vectors_pp.dat', status='unknown')
              
              n_aux=1
              do i = 1, Nt
                if (dble(EigenVa(i)) > Zero .and. dble(EigenVa(i)) < One+1.d-2) then
                  
                  do iw = 1,Nf
                    do ix2 = 1,Nx
                      do iy2= 1,Ny
                        idx = ((ix2-1)*Ny+iy2-1)*Nf + iw
                        
                        write(133, '(1i5, 1f15.8, 3i5, 2f15.8)') i, dble(EigenVa(i)), ix2, iy2, (iw-Nf/2-1), dble(EigenVR(idx,i)), aimag(EigenVR(idx,i))
                      end do
                    end do
                  end do
                  
                  n_aux=n_aux+1
                end if
                
                if (n_aux>20) exit
              end do
              
              close(133)  
              
              ! call ZGEMM('N', 'N', Nt, Nt, Nt, One_c, EigenVR, Nt, dummy3D_1, Nt, One_c, EigenVL, Nt)

!! Obsolete in version 1.1
              
!               do i = 1,Nt
!                 map_i =index_fermionic(i)
!                 do j = 1,Nt            
!                   map_j = index_fermionic(j)
!                   
!                   Susc_s = Susc_s + dummy3D_1(i,j)
!                   Susc_p = Susc_p + gsym_p(map_i%ix)*gsym_p(map_j%ix)*dummy3D_1(i,j)
!                   Susc_dxy = Susc_dxy + gsym_dxy(map_i%ix,map_i%iy)*gsym_dxy(map_j%ix,map_j%iy)*dummy3D_1(i,j)
!                   Susc_dx2y2 = Susc_dx2y2 + gsym_dx2y2(map_i%ix,map_i%iy)*gsym_dx2y2(map_j%ix,map_j%iy)*dummy3D_1(i,j)
!                   
!                 end do
!               end do            
!               
!               Susc_s = Susc_s/beta/beta/Nc
!               Susc_p = Susc_p/beta/beta/Nc
!               Susc_dxy = Susc_dxy/beta/beta/Nc              
!               Susc_dx2y2 = Susc_dx2y2/beta/beta/Nc
              
              if (allocated(dummy3D_1)) deallocate(dummy3D_1)
              
              !  write(FLE2, '(a)') 'EigenVal_ite_pp.dat'
              open(unit=130, file=FLE2, status='unknown', position='append')
              
              n_aux=1
              do i = 1, Nt
                if (dble(EigenVa(i)) > Zero) then
                  write(130, '(1i5, 2f15.8)') ite, dble(EigenVa(i)), aimag(EigenVa(i))
                  !, dble(Susc_s), dble(Susc_dx2y2),dble(Susc_dxy), dble(Susc_p)
                  
                  if (n_aux==1) eign_pp=dble(EigenVa(i))
                  n_aux=n_aux+1
                end if
                if (n_aux>10) exit
              end do
              
              close(130) 
              
            end if  !ichannel >= 3  
            
            if ((map_q%iy == Nx_IBZ).AND.(iChannel == 2)) then 
              
              !    write(FLE2, '(a)') 'EigenVal_ite_m.dat'             
              
              open(unit=131, file=FLE2, status='unknown', position='append')
              
              n_aux=1
              do i = 1, Nt
                if (dble(EigenVa(i)) > Zero) then
                  write(131, '(1i5, 6f15.8)') ite, dble(EigenVa(i)), aimag(EigenVa(i))
                  
                  if (n_aux==1) eign_m=dble(EigenVa(i))
                  n_aux=n_aux+1
                end if
                
                if (n_aux>10) exit
              end do
              close(131) 
              
              if (.NOT. allocated(dummy3D_1)) allocate(dummy3D_1(Nt, Nt))
              
              ! call ZGEMM('C', 'N', Nt, Nt, Nt, One_c, EigenVL, Nt, EigenVR, Nt, Zero_c, dummy3D_1, Nt)
              
              do i = 1,Nt
                do j = 1,Nt
                  dummy3D_1(i,j)=-dummy1D(i)*(F_m(i,j,k))*dummy1D(j)
                end do
              end do
              
              open(unit=123, file='Susc_m.dat', status='unknown')
              
              do ix2 = 1,Nx
                do iy2= 1,Ny
                  do iw2 = 1,Nf
                    idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
                    do ix = 1,Nx
                      do iy= 1,Ny
                        do iw = 1,Nf
                          idx = ((ix-1)*Ny+iy-1)*Nf + iw       
                          
                          write(123, '( 6i5, 2f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(dummy3D_1(idx2,idx)), aimag(dummy3D_1(idx2,idx))
                        end do
                      end do
                    end do
                  end do
                end do
              end do             
              
              close(123)  
              
              call output_2k(F_d, F_m, k, 'F_ph_Q.dat', 127)
!               open(unit=127, file='F_ph_Q.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           write(127, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(F_d(idx2,idx,k)), aimag(F_d(idx2,idx,k)), dble(F_m(idx2,idx,k)), aimag(F_m(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(127) 
              
              call output_2k(F_s, F_t, k, 'F_pp_Q.dat', 128)
!               open(unit=128, file='F_pp_Q.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(128, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(F_s(idx2,idx,k)), aimag(F_s(idx2,idx,k)), dble(F_t(idx2,idx,k)), aimag(F_t(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(128) 
              
              call output_2k(G_d, G_m, k, 'Gamma_ph_Q.dat', 129)
!               open(unit=129, file='Gamma_ph_Q.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(129, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(G_d(idx2,idx,k)), aimag(G_d(idx2,idx,k)), dble(G_m(idx2,idx,k)), aimag(G_m(idx2,idx,k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!               
!               close(129) 
              
              call output_2k(G_s, G_t, k, 'Gamma_pp_Q.dat', 122)
              ! Technically, this might break some scripts for post-processing, but
              !call output_2k(0.5*(G_t-G_s), G_s, 'Gamma_pp_Q.dat', 122)
              ! should work. It might however, cause a stack overflow
!               open(unit=122, file='Gamma_pp_Q.dat', status='unknown')
!               
!               do ix2 = 1,Nx
!                 do iy2= 1,Ny
!                   do iw2 = 1,Nf
!                     idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
!                     do ix = 1,Nx
!                       do iy= 1,Ny
!                         do iw = 1,Nf
!                           idx = ((ix-1)*Ny+iy-1)*Nf + iw       
!                           
!                           write(122, '(6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1,&
!                           dble(0.5d0*(G_t(idx2, idx, k)-G_s(idx2, idx, k))), aimag(0.5d0*(G_t(idx2, idx, k)-G_s(idx2, idx, k))), dble(G_s(idx2, idx, k)), aimag(G_s(idx2, idx, k))
!                         end do
!                       end do
!                     end do
!                   end do
!                 end do
!               end do             
!              
!               close(122)  
              
            end if !ichannel == 2
            
            if(Converged) then
              
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
                if (dble(EigenVa(i)) > Zero .and. dble(EigenVa(i)) < One+1.d-2) then
                  !select case (ichannel)
                    !case (1,2)
                      write(132, '(2i5, 2f15.8)') map_q%ix, map_q%iy, dble(EigenVa(i)), aimag(EigenVa(i))
                      
                    !case (3)
! ! Obsolete in version 1.1
!                       Ps = Zero
!                       Pp = Zero
!                       Pdxy = Zero
!                       Pdx2y2 = Zero
!                       
!                       do ix2 = 1,Nx
!                         do iy2= 1,Ny
!                           do iw = Nf/2+1,Nf
!                             idx = ((ix2-1)*Ny+iy2-1)*Nf + iw
!                             Ps = Ps + EigenVR(idx,i)
!                             Pp = Pp + EigenVR(idx,i)*gsym_p(ix2)
!                             Pdxy = Pdxy + EigenVR(idx,i)*gsym_dxy(ix2,iy2)
!                             Pdx2y2 = Pdx2y2 + EigenVR(idx,i)*gsym_dx2y2(ix2,iy2)
!                           end do
!                         end do
!                       end do
                      !write(132, '(2i5, 2f15.8)') map_q%ix, map_q%iy, dble(EigenVa(i)), aimag(EigenVa(i))
                      !, abs(Ps)/beta/Nc, abs(Pdx2y2)/beta/Nc,abs(Pdxy)/beta/Nc,abs(Pp)/beta/Nc
                  !end select
                  
                end if
              end do
              
              close(132)          
              
            end if ! Converged or not
          end if ! Channel and qx,qy
        end do ! iChannel
      end if !qw = 0
    end do !Nb
  end if ! Eigen  
    
    
    if (id == master) then 
      
      open(unit=102, file='Error_ite.dat', status='unknown', position='append') 
      
      write(102, '(1i5, 4f12.6)') ite, relative_err, mu, eign_m,  eign_pp
      
      close(102) 
      
    end if 
    
  end subroutine solve_eigen_equation
  
  ! -----------------------------------------------------------------------------------------------------------
end module Parquet_Phys
