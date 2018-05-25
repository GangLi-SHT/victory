module parquet_selfenergy

  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util
  use parquet_kernel
  use parquet_equation

contains
  !----------------------------------------------------------------------
  subroutine self_energy(ite, Grt, Converged)
    !
    ! Purpose
    ! ========
    !  calculate self-energy from Schwinger-Dyson equation. 
    !  The calculated self-energy is then distributed to each node. 
    !
    implicit none
    integer, intent(in)     :: ite
    complex(dp), intent(in) :: Grt(Nx, Ny, Nf)
    logical, intent(out)    :: Converged

    ! ... local vars ...
    integer       :: idx, iTau
    integer       :: i, j, k, i1, j1, jx, jy, jw
    integer       :: ichannel
    real(dp)      :: relative_err, FD1, FD2, dkx, dky, kx, ky
    complex(dp)   :: dummy, dummy1, dummy2, dummy3, dummy4, Sigma_V(Nt)
    complex(dp), allocatable :: SigmaOld(:), dummy3d_1(:,:,:), coutdata(:), Sigma_Reduced(:)
    type(Indxmap) :: ComIdx1, ComIdx2, ComIdx3, ComIdx4, ComIdx5, ComIdx_j  ! combined index for momentum and frequency
    type(Indxmap) :: map_i, map_j, map_k
    character(len=30) :: FLE, str1, str2
    character(len=10) :: Mtype

    if (.NOT. allocated(SigmaOld)) allocate(SigmaOld(Nt))
    SigmaOld = Sigma
    Sigma = Zero

    ! --- determine the rest self-energy diagrams in the 2nd parquet approximation ---

    !  [1] Hartree term, which becomes nonzero at away-halfing case. At
    !  half-filling this term is U/2 and is absorbed into the chemical
    !  potential, thus, the chemical potential for half-filled case in this code is
    !  defined as zero. 
    do i = 1, Nx
       do j = 1, Ny
          do k = 1, Nf
             idx = (Ny*(i-1)+j-1)*Nf + k
             Sigma(idx) = -xU*(Grt(1, 1, Nf) + Half)
          end do
       end do
    end do
    
    ! [2] second order diagram 
    
    if (.NOT. allocated(dummy3d_1)) allocate(dummy3d_1(Nx, Ny, Nf))
    MType = 'Fermionic'
    dummy3D_1 = Zero
    do iTau = 1, Nf
       do i = 1, Nx
          if (i == 1) then
             i1 = 1
          else
             i1 = Nx-i+2
          end if
          do j = 1, Ny
             if (j == 1) then
                j1 = 1
             else
                j1 = Ny-j+2
             end if
             dummy3D_1(i, j, iTau) = Grt(i, j, iTau)**2*Grt(i1, j1, Nf-iTau+1)*xU*xU
          end do
       end do
       
       ! --- now change it to k-t space by using FFT ---
       call fftb2d(Nx, Ny, dummy3D_1(1:Nx, 1:Ny, iTau), C_wave_x, C_wave_y)    ! 2nd order Sigma in k-t space!       
    end do
    
    ! 2nd order Sigma in k-w space
    if (.NOT. allocated(coutdata)) allocate(coutdata(Nf/2))
    do i = 1, Nx
       do j = 1, Ny
          call FDfit(Nf, dble(dummy3D_1(i, j, 1:Nf)), beta/(Nf-1), FD1, FD2)
          call nfourier(MType, Nf-1, Nf/2, FD1, FD2, dble(dummy3D_1(i, j, 1:Nf)), coutdata)
          do k = 1, Nf/2
             idx = (Ny*(i-1)+j-1)*Nf + k
             Sigma(idx) = Sigma(idx) + conjg(coutdata(Nf/2+1-k))
             Sigma(idx+Nf-2*k+1) = Sigma(idx+Nf+1-2*k) + coutdata(Nf/2+1-k)
          end do
       end do
    end do
    
    if (id == master) then
       ! --- Output the self-energy function from 2nd parquet approximation ---
       write(str1, '(I0.3)') ite
       FLE = 'Sigma2nd-'//trim(str1)//'.dat'
       open(unit=1, file=FLE, status='unknown')
       dkx = Two*Pi/dble(Nx)
       dky = Two*Pi/dble(Ny)
       do i = 1, Nx
          kx = dkx*(i-1)
          do j = 1, Ny
             ky = dky*(j-1)
             do k = 1, Nf
                idx = (Ny*(i-1)+j-1)*Nf + k
                write(1, '(7f20.12)') kx, ky, Pi/beta*(Two*(k-Nf/2-1)+One), Sigma(idx)
             end do
             write(1, *)
          end do
       end do
       close(1)
    end if
    if (allocated(dummy3d_1)) deallocate(dummy3d_1)
    if (allocated(coutdata))  deallocate(coutdata)
      
    ! the rest of the self-energy
    Sigma_V = Zero
    do ichannel = 1, 2
       ! denisty and magnetic channels are sufficient for constructing self-energy. 
       select case (ichannel)
       case (1)
          mat = F_d - xU
       case (2)
          mat = F_m + xU
       end select
       
       ! in the following, the self-energy is only determined inside the box using the full-vertex. 
       do i = 1, Nt
          dummy = Zero

          map_i = index_fermionic(i)
          call index_operation(map_i, map_i, MinusF, ComIdx1)
          ! variable k'

          do j = 1, Nt

          map_j = index_fermionic(j)
           dummy3 = Gkw(j)
           
             ! variable q             
             do k = 1, Nb
            map_k = index_bosonic(id*Nb+k)

            call index_operation(map_i, map_k, FaddB, ComIdx2) ! k+q
                if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
                   dummy1 = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) )
                else
                   dummy1 = Gkw(list_index(ComIdx2, Fermionic))
                end if
                call index_operation(map_j, map_k, FaddB, ComIdx3) ! k'+q
                if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                   dummy2 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) )
                else
                   dummy2 = Gkw(list_index(ComIdx3, Fermionic))
                end if
                
                   dummy = dummy + dummy3*dummy1*dummy2*(mat(i, j, k))

                if (map_k%iw > 1) then
                   ! [2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                   
                   ! -U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(k'+q)^* x G(k-q) x G(k')^* ] 
                   call index_operation(ComIdx1, map_k, FaddB, ComIdx4)  ! -k+q
                   call index_operation(ComIdx4, ComIdx4, MinusF, ComIdx5) ! k-q
                   if (ComIdx5%iw > Nf .or. ComIdx5%iw < 1) then
                      dummy1 = One/(xi*Pi/beta*(Two*(ComIdx5%iw-Nf/2-1) + One) + mu - Ek(ComIdx5%ix, ComIdx5%iy) )
                   else
                      dummy1 = Gkw(list_index(ComIdx5, Fermionic))
                   end if
                      dummy = dummy + conjg(dummy3)*dummy1*conjg(dummy2)*conjg(mat(list_index(ComIdx1, Fermionic), j, k))
                end if
             end do
          end do
                   select case (ichannel)
                   case(1)
          Sigma_V(i) = Sigma_V(i) - xU/beta/beta/Nc/Nc/2.0*dummy
                    case(2)
          Sigma_V(i) = Sigma_V(i) + xU/beta/beta/Nc/Nc/2.0*dummy                    
                    
                    end select         
       end do
              
    end do

    ! now we determine the self-energy solely from the auxiliary function in an enlarged parameter space
    do i = 1, Nt
       dummy = 0.d0
       map_i=index_fermionic(i)
       call index_operation(map_i, map_i, MinusF, ComIdx1)
       do jx = 1, Nx
          do jy = 1, Ny
             do jw = -Nf+1, 2*Nf
              
                ! Gkw(k')
                if (jw > Nf .or. jw < 1) then
                  ComIdx_j = indxmap(jx, jy, jw)
 
                dummy1 = One/(xi*Pi/beta*(Two*(jw-Nf/2-1) + One) + mu - Ek(jx, jy) )
                
                do k = 1, Nb
                    map_k= index_bosonic(id*Nb+k)
                    ! Gkw(k'+q)
                   call index_operation(ComIdx_j, map_k, FaddB, ComIdx2) ! k'+q

                   if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
                      dummy2 = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) )
                   else
                      dummy2 = Gkw(list_index(ComIdx2, Fermionic))
                   end if

                   ! Gkw(k+q)
                   call index_operation(map_i, map_k, FaddB, ComIdx3) ! k+q

                   if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                      dummy3 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) )
                   else
                      dummy3 = Gkw(list_index(ComIdx3, Fermionic))
                   end if
                   
                   dummy = dummy - dummy1*dummy2*dummy3*( &
                   kernel('d', map_i, ComIdx_j,map_k) - 3.0d0*kernel('m', map_i, ComIdx_j, map_k))
                                     
                
                   ! [2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                   
                   ! -U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(k'+q)^* x G(k-q) x G(k')^* ] 
                                      
                    call index_operation(ComIdx1, map_k, FaddB, ComIdx3)  ! -k+q
                    call index_operation(ComIdx3, ComIdx3, MinusF, ComIdx4) ! k-q               
                   
                                      
                   if ( map_k%iw > 1) then
                      ! time-reversal part 

                      if (ComIdx4%iw > Nf .or. ComIdx4%iw < 1) then
                         dummy4 = One/(xi*Pi/beta*(Two*(ComIdx4%iw-Nf/2-1) + One) + mu - Ek(ComIdx4%ix, ComIdx4%iy) )
                      else
                         dummy4 = Gkw(list_index(ComIdx4, Fermionic))
                      end if
                      
                      
                      dummy = dummy - conjg(dummy1)*conjg(dummy2)*dummy4*conjg(&
                      kernel('d', ComIdx1, ComIdx_j,map_k) - 3.0d0*kernel('m', ComIdx1, ComIdx_j,map_k))
           
                 end if
                 
                 ! The contributions of s and t kernels have shifted variables: k+k'+q -> q, k+q -> q-k', k'+q -> q-k   
  
                 ! Gkw(q-k')
                    call index_operation(ComIdx_j, ComIdx_j, MinusF, ComIdx5) ! -k'
                    call index_operation(ComIdx5, index_bosonic(id*Nb+k), FaddB, ComIdx2) ! -k'+q
                   if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
                      dummy2 = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) )
                   else
                      dummy2 = Gkw(list_index(ComIdx2, Fermionic))
                   end if

                   ! Gkw(q-k)
                   if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then  ! -k+q
                      dummy4 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) )
                   else
                      dummy4 = Gkw(list_index(ComIdx3, Fermionic))
                   end if
                   
                   
                    dummy = dummy - dummy1*dummy2*dummy4*( &
                    kernel('s', map_i, ComIdx_j, map_k) + kernel('t', map_i, ComIdx_j,map_k))
                 
                  if ( index_bosonic(id*Nb+k)%iw > 1) then
                                 
                         ! [2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                   ! -U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(q-k')^* x G(-k-q) x G(k')^* ] =
                   !  -U x T^2/2/N^2 sum_k' sum_q conjg {[Fd(-k,k';q) - Fm(-k,k';q)] x G(q-k') x G(k+q) x G(k') }
                   
                    dummy = dummy - conjg(dummy1*dummy2*dummy3*( &
                    kernel('s',  ComIdx1, ComIdx_j, map_k) + kernel('t', ComIdx1, ComIdx_j,map_k )))                
                 
                  
                   end if
                end do
                end if
             end do
          end do
       end do
       Sigma_V(i) = Sigma_V(i) +  0.5d0*xU/beta/beta/Nc/Nc*dummy
    end do
    
    IF (.NOT. ALLOCATED(Sigma_Reduced)) ALLOCATE(Sigma_Reduced(Nt))

    call MPI_AllReduce(Sigma_V, Sigma_Reduced, Nt, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
    Sigma = Sigma + Sigma_Reduced
    
    Deallocate(Sigma_Reduced)

    ! output the Green's function in k-w space
    if (id == master) then 
       write(str1, '(I0.3)') ite
       FLE = 'Sigmakw-'//trim(str1)//'.dat'
       open(unit=1, file=FLE, status='unknown')
       do i = 1, Nx
          do j = 1, Ny
             do k = 1, Nf
                idx = (Ny*(i-1)+j-1)*Nf + k
                write(1, '(2i4, 3f20.12)') i, j, Pi/beta*(Two*(k-Nf/2)-One), Sigma(idx)
             end do
             write(1, *)
          end do
       end do
       close(1)

       ! --- checking convergence ---
       relative_err = Zero
       do i = 1, Nt
          relative_err = relative_err + abs(One-SigmaOld(i)/Sigma(i))/Nt
       end do
       if (relative_err < 1.d-3) then
          Converged = .True.
       else
          Converged = .False.
       end if
       write(*, '(a, f12.6)') 'relative error in self-energy is:', relative_err
    end if

    call MPI_BCAST(Converged, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, rc)

    if (Converged) then
       write(str1, '(I0.3)') id
       write(str2, '(I0.3)') ite
       FLE = 'F-'//trim(str1)//'-'//trim(str2)//'.dat'
       open(unit=1, file=FLE, status='unknown')
       do i = 1, Nt
          do j =1, Nt
             do k = 1, Nb
                write(1, '(3i5, 8f12.6)') i, j, id*Nb+k, F_d(i, j, k), F_m(i, j, k), F_s(i, j, k), F_t(i, j, k)
             end do
          end do
       end do
       close(1)
    end if

    if (id == master) call system_mem_usage
  end subroutine self_energy

end module parquet_selfenergy
