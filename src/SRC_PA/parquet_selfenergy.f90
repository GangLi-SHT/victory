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
    !  calculate self-energy from Schwinger-Dyson equation. The self-energy is only calculated at master node, this reduces the
    ! communication between different nodes to the minimum. The calculated self-energy is then distributed to each node afterwards. 
    !
    implicit none
    integer, intent(in)     :: ite
    complex(dp), intent(in) :: Grt(Nx, Ny, Nf)
    logical, intent(out)    :: Converged
    
    ! ... local vars ...
    integer       :: idx, iTau
    integer       :: i, j, k, i1, j1, jx, jy, jw, is, ix, iy, iw, icx, icy, jx2, jy2, j3
    !integer       :: ichannel
    real(dp)      :: FD1, FD2, dkx, dky, kx, ky, w, wc
    complex(dp)   :: dummy, dummy1, dummy2, dummy3, dummy4, Sigma_V(Nred*2)
    complex(dp), allocatable :: SigmaOld(:), dummy3d_1(:,:,:), coutdata(:), Sigma_Reduced(:)
    complex(dp), allocatable :: dummy2d_1(:,:,:,:,:), dummy2d_2(:,:,:,:,:)
    type(Indxmap) :: ComIdx1, ComIdx2, ComIdx3, ComIdx4, ComIdx5, ComIdx6   ! combined index for momentum and frequency
    type(Indxmap) :: map_i, map_j, map_k, map_i1, map_j1, map_k1
    character(len=30) :: FLE, str1, str2
    character(len=10) :: Mtype
    logical :: symm_list(Ns)
    
    if (.NOT. allocated(SigmaOld)) allocate(SigmaOld(Nt))
    SigmaOld = Sigma
    Sigma = Zero
    
    ! --- determine the rest self-energy diagrams in the 2nd parquet approximation ---
    
    !  [1] Hartree term, which becomes nonzero at away-halfing case. At
    !  half-filling this term is U/2 and is absorbed into the chemical
    !  potential, thus, the chemical for half-filled case in this code is
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
      call fftb2d(Nx, Ny, dummy3D_1(1:Nx, 1:Ny, iTau), C_wave_x, C_wave_y)          ! 2nd Sigma in k-t space!       
    end do
    
    ! 2nd Sigma in k-w space
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
    
    ! Auxiliary calculation of G(k')G(k'+q) for the self-energy with coarse graining
    
    if (.NOT. allocated(dummy2d_1)) allocate(dummy2d_1(Ny, Nx, (2*f_range+1)*Nf,Ns,Nb))
    if (.NOT. allocated(dummy2d_2)) allocate(dummy2d_2(Ny, Nx, (2*f_range+1)*Nf,Ns,Nb))
    
    dummy2d_1=Zero
    dummy2d_2=Zero
    
    ! variable q             
    do k = 1, Nb !q
      map_k= index_bosonic_IBZ(id*Nb+k)
      call list_symmetries(lattice_type,map_k,symm_list)
      
      !symmetries    
      do is= 1, Ns ! the number of operations is lattice dependent
        
        if (symm_list(is)) then
          ! is=1
          call symmetry_operation(lattice_type,is,map_k, map_k1) !S(q)
          
          ! variable k'
          do jw = -f_range*Nf+1, (f_range+1)*Nf !k'
            
            do jx = 1, Nx
              do jy = 1, Ny              
                map_j = indxmap(jx, jy, jw) !k' do j = 1, Nt
                map_i = indxmap(ix, iy, iw)
                
                j=list_index_F(map_j)
                
                call symmetry_operation_inv(lattice_type,is,map_i,map_j, map_i1, map_j1) !S^-1(k')
                
                j1=list_index_F(map_j1)
                
                call index_FaddB(map_j, map_k1,  ComIdx3) ! k'+S(q)
                
                j3=list_index_F(ComIdx3)
                
                icx = ComIdx3%ix
                icy = ComIdx3%iy
                
                
                dummy2 = Zero
                dummy3 = Zero
                
                w = Pi/beta*(Two*(jw-Nf/2-1) + One)
                wc = Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One)
                
                if (jw > Nf .or. jw < 1) then
                  
                  if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        dummy2 = dummy2 + One/(xi*wc + mu - Ek_grain(icx,jx2, icy, jy2) -&
                        Sigma_H(icx,icy))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2)-Sigma_H(jx,jy))
                      end do
                    end do
                  else
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        dummy2 = dummy2 + One/(xi*wc + mu - Ek_grain(icx,jx2, icy, jy2) - &
                        SigmaOld(j3))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2)-Sigma_H(jx,jy)) 
                      end do
                    end do
                    !Gkw(list_index_F(ComIdx2))
                  end if 
                  
                  
                ! Gkw(q-k')
                  call index_minusF(map_j,  ComIdx5) ! -k' 
                  call index_FaddB(ComIdx5, map_k1, ComIdx2) ! -k'+S(q)
                  
                  icx = ComIdx2%ix
                  icy = ComIdx2%iy
                  j3=list_index_F(ComIdx2)
                  
                  wc = Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One)
                  
                  if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        
                        dummy3 = dummy3 + One/(xi*wc + mu - Ek_grain(icx,Ngrain-jx2+1, icy,Ngrain- jy2+1) - &
                        Sigma_H(icx,icy))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2)-Sigma_H(jx,jy))
                      end do
                    end do
                  else
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        dummy3 = dummy3 + One/(xi*wc + mu - Ek_grain(icx,Ngrain-jx2+1, icy, Ngrain-jy2+1) - &
                        SigmaOld(j3))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2)-Sigma_H(jx,jy)) 
                      end do
                    end do
                    !Gkw(list_index_F(ComIdx2))
                  end if                           
                  
                  
                end if
                
                if (jw <= Nf .and. jw >= 1) then
                  
                  if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        dummy2 = dummy2 + One/(xi*wc + mu - Ek_grain(icx,jx2, icy, jy2) - &
                        Sigma_H(icx, icy))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2) - SigmaOld(j))
                      end do
                    end do
                  else
                    do jx2 = 1,Ngrain
                      do jy2 = 1,Ngrain
                        dummy2 = dummy2 + One/(xi*wc + mu - Ek_grain(icx,jx2, icy, jy2) - &
                        SigmaOld(j3))/(xi*w + mu - Ek_grain(jx,jx2, jy, jy2) - SigmaOld(j)) 
                      end do
                    end do
                   ! Gkw(list_index_F(ComIdx3))
                  end if 
                  
                end if   
               
                dummy2d_1(jy,jx,jw+f_range*Nf,is,k) = dummy2/Ngrain/Ngrain
                dummy2d_2(jy,jx,jw+f_range*Nf,is,k) = dummy3/Ngrain/Ngrain
                
              end do !k'
            end do !k'
            
          end do !k'
        end if !symm_list
      end do ! symmetries
      
    end do ! q
    
    ! the rest of the self-energy
    Sigma_V = Zero
    ! denisty and magnetic channels are sufficient for constructing self-energy. 
    mat = F_d - F_m - Two*xU
    
    ! in the following, the self-energy is only determined inside the box by subtracting the full-vertex with the auxiliary funcition. 
    do ix = 1, Nx_IBZ  ! variable k
      do iy = 1,ix
        do iw = 1, Nf
          
          dummy = Zero
          
          map_i%ix = ix
          map_i%iy = iy
          map_i%iw = iw
          
          i = (((ix-1)*ix)/2+iy-1)*Nf + iw !list_index(map_i,Fermionic)
          
          call index_minusF(map_i,  ComIdx6)
          
          ! variable q             
          do k = 1, Nb !q
            map_k= index_bosonic_IBZ(id*Nb+k)
            call list_symmetries(lattice_type,map_k,symm_list)
            
            !symmetries    
            do is= 1, Ns ! the number of operations is lattice dependent
              
              if (symm_list(is)) then
                !is=1
                call symmetry_operation(lattice_type,is,map_k, map_k1) !S(q)
                
                ! variable k'
                do j = 1, Nt
                  map_j=index_fermionic(j)
                  !dummy3 = Gkw(j) 
                  
                  call symmetry_operation_inv(lattice_type,is,map_i,map_j, map_i1, map_j1) !S^-1(k),S^-1(k')
                  
                  i1=list_index_F(map_i1)
                  j1=list_index_F(map_j1)
                  
                  call index_minusF(map_i1,  ComIdx1)
                  
                  call index_FaddB(map_i, map_k1, ComIdx2) ! k+S(q)
                  
                  if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
                    dummy1 = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy)- Sigma_H(ComIdx2%ix, ComIdx2%iy) )
                  else
                    dummy1 = One/(xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) - SigmaOld(list_index_F(ComIdx2))) 
                  end if
                  
                  !call index_FaddB(map_j, map_k1,  ComIdx3) ! k'+S(q)
                  
                  jx = map_j%ix
                  jy = map_j%iy 
                  jw = map_j%iw
                  
                  dummy2 = dummy2d_1(jy,jx,jw+f_range*Nf,is,k)
                  dummy = dummy + dummy1*dummy2*(mat(i1, j1, k))
                  
                  if (map_k%iw > 1) then
                    ! [2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                    
                    ! -U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(k'+q)^* x G(k-q) x G(k')^* ] 
                    call index_FaddB(ComIdx6, map_k1,ComIdx4)  ! -k+S(q)
                    call index_minusF(ComIdx4, ComIdx5) ! k-S(q)
                   
                    if (ComIdx5%iw > Nf .or. ComIdx5%iw < 1) then
                      dummy1 = One/(xi*Pi/beta*(Two*(ComIdx5%iw-Nf/2-1) + One) + mu - Ek(ComIdx5%ix, ComIdx5%iy) - Sigma_H(ComIdx5%ix, ComIdx5%iy) )
                    else
                      dummy1 = One/(xi*Pi/beta*(Two*(ComIdx5%iw-Nf/2-1) + One) + mu - Ek(ComIdx5%ix, ComIdx5%iy) - SigmaOld(list_index_F(ComIdx5)))!Gkw(list_index_F(ComIdx5))
                    end if
                    
                    dummy = dummy + dummy1*conjg(dummy2)*conjg(mat(list_index_F(ComIdx1), j1, k))
                  end if
                end do ! k'
              end if !symm_list
            end do !symmetries
          end do !q
          
          Sigma_V(i) = Sigma_V(i) - xU/beta/beta/Nc/Nc/2.0*dummy
          
        end do !k
      end do !k
    end do !k
    
    !end do !channels
    !now determine the self-energy solely from the auxiliary function in an enlarged parameter space
    !do i = 1, Nt !k
    
    
    do ix = 1, Nx_IBZ  ! variable k
      do iy = 1,ix
        do iw = 1, Nf
          dummy = 0.d0       
          map_i%ix = ix
          map_i%iy = iy
          map_i%iw = iw
          
          i = (((ix-1)*ix)/2+iy-1)*Nf + iw !list_index(map_i,Fermionic)
          
          call index_minusF(map_i, ComIdx6)
          
          do k = 1, Nb !q
            map_k= index_bosonic_IBZ(id*Nb+k)
            call list_symmetries(lattice_type,map_k,symm_list)
            
            !symmetries    
            do is= 1, Ns ! the number of operations is of course lattice dependent
              
              if (symm_list(is)) then
                
                call symmetry_operation(lattice_type,is,map_k, map_k1) !S(q)      
                
                do jw = -f_range*Nf+1, (f_range+1)*Nf !k'
                  if (jw > Nf .or. jw < 1) then
                    do jx = 1, Nx
                      do jy = 1, Ny              
                        map_j = indxmap(jx, jy, jw) !k'
                        
                        !Gkw(k')
                        !dummy1 = One/(xi*Pi/beta*(Two*(map_j%iw-Nf/2-1) + One) + mu - Ek(map_j%ix, map_j%iy) )
                        
                        call index_minusF(map_j,  ComIdx5) ! -k'
                        
                        call symmetry_operation_inv(lattice_type,is,map_i,map_j, map_i1, map_j1) !!S^-1(k),S^-1(k')
                        
                        i1=list_index_F(map_i1)
                        j1=list_index_F(map_j1)
                        
                        call index_minusF(map_i1, ComIdx1)               
                        !Gkw(k'+q)
                        !call index_FaddB(map_j, map_k1, ComIdx2) ! k'+S(q)
                        
                        dummy2 = dummy2d_1(jy,jx,jw+f_range*Nf,is,k)
                        
                        !Gkw(k+q)
                        call index_FaddB(map_i, map_k1, ComIdx3) ! k+S(q)
                        
                        if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then
                          dummy3 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy)- Sigma_H(ComIdx3%ix, ComIdx3%iy))
                        else
                          dummy3 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) - SigmaOld(list_index_F(ComIdx3)))!Gkw(list_index_F(ComIdx3))
                        end if
                        
                        dummy = dummy - dummy2*dummy3*( &
                        kernel('d', map_i1, map_j1,map_k) - 3.0d0*kernel('m', map_i1, map_j1, map_k))
                        
                        ![2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                        !-U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(k'+q)^* x G(k-q) x G(k')^* ] 
                        
                        call index_FaddB(ComIdx6, map_k1,  ComIdx3)  ! -k+S(q)
                        call index_minusF(ComIdx3, ComIdx4) ! k-S(q)               
                        
                        if ( map_k%iw > 1) then
                          !time-reversal part 
                          if (ComIdx4%iw > Nf .or. ComIdx4%iw < 1) then
                            dummy4 = One/(xi*Pi/beta*(Two*(ComIdx4%iw-Nf/2-1) + One) + mu - Ek(ComIdx4%ix, ComIdx4%iy) - Sigma_H(ComIdx4%ix, ComIdx4%iy))
                          else
                            dummy4 = One/(xi*Pi/beta*(Two*(ComIdx4%iw-Nf/2-1) + One) + mu - Ek(ComIdx4%ix, ComIdx4%iy) - SigmaOld(list_index_F(ComIdx4)))! Gkw(list_index_F(ComIdx4))
                          end if
                          
                          dummy = dummy - conjg(dummy2)*dummy4*conjg(&
                          kernel('d', ComIdx1, map_j1,map_k) - 3.0d0*kernel('m', ComIdx1, map_j1,map_k))
                        end if
                        
                        !The contributions of s and t kernels have shifted variables: k+k'+q -> q, k+q -> q-k', k'+q -> q-k   
                        !Gkw(q-k')
                        
                        !call index_FaddB(ComIdx5, map_k1, ComIdx2) ! -k'+S(q)
                        dummy2 = dummy2d_2(jy,jx,jw+f_range*Nf,is,k)
                        
                        ! Gkw(q-k)
                        if (ComIdx3%iw > Nf .or. ComIdx3%iw < 1) then  ! -k+S(q)
                          dummy4 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) - Sigma_H(ComIdx3%ix, ComIdx3%iy))
                        else
                          dummy4 = One/(xi*Pi/beta*(Two*(ComIdx3%iw-Nf/2-1) + One) + mu - Ek(ComIdx3%ix, ComIdx3%iy) - SigmaOld(list_index_F(ComIdx3)))!Gkw(list_index_F(ComIdx3))
                        end if
                        
                        dummy = dummy - dummy2*dummy4*( &
                        kernel('s', map_i1, map_j1, map_k) + kernel('t', map_i1, map_j1,map_k))
                        
                        if ( map_k%iw > 1) then
                          ! [2]: time-reversal symmetry conjg[ F(k,k';q) ] = F(-k,-k';-q) 
                          ! -U x T^2/2/N^2 sum_k' sum_q [Fd(-k,k';q) - Fm(-k,k';q)]^* x [G(q-k')^* x G(-k-q) x G(k')^* ] =
                          ! -U x T^2/2/N^2 sum_k' sum_q conjg {[Fd(-k,k';q) - Fm(-k,k';q)] x G(q-k') x G(k+q) x G(k') }
                        
                          dummy = dummy - conjg(dummy2*dummy3*( &
                          kernel('s',  ComIdx1, map_j1, map_k) + kernel('t', ComIdx1, map_j1,map_k )))
                          
                        end if
                      end do !k'
                    end do !k'
                  end if
                end do !k'
                
              end if !symm_list
            end do ! symmetries
            
          end do ! q
          Sigma_V(i) = Sigma_V(i) +  0.5d0*xU/beta/beta/Nc/Nc*dummy
        end do !k
      end do
    end do
    
    if (allocated(dummy2d_1)) deallocate(dummy2d_1)
    if (allocated(dummy2d_2)) deallocate(dummy2d_2)
    
    if (.NOT. ALLOCATED(Sigma_Reduced)) allocate(Sigma_Reduced(Nred*2))
    
    call MPI_AllReduce(Sigma_V, Sigma_Reduced, Nred*2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
    
    do ix = 1, Nx_IBZ  ! variable k
      do iy = 1,ix
        do iw = 1, Nf
          
          map_i%ix = ix
          map_i%iy = iy
          map_i%iw = iw
          
          i = (((ix-1)*ix)/2+iy-1)*Nf + iw
          
          call list_symmetries(lattice_type,map_i,symm_list)
          !symmetries    
          do is= 1, Ns ! the number of operations is of course lattice dependent          
            if (symm_list(is)) then
              
              call symmetry_operation(lattice_type,is,map_i, map_i1) !S(k)           
              
              i1=list_index_F(map_i1)
              Sigma(i1) = Sigma(i1) + Sigma_Reduced(i)
              
            end if
          end do
        end do
      end do
    end do
    
    deallocate(Sigma_Reduced)
    
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
      
      if(.NOT.No_sigma_update) then
        do i = 1, Nt
          relative_err = relative_err + abs(One-SigmaOld(i)/Sigma(i))/Nt
        end do
      else
        do i = 1, Nt
          relative_err = relative_err + abs(One-Sigma_compare(i)/Sigma(i))/Nt
        end do
      end if
      
      if (relative_err < f_sigma) then
        Converged = .True.
      else
        Converged = .False.
      end if
      write(*, '(a, f12.6)') 'relative error in self-energy is:', relative_err
    end if
    
    call MPI_BCAST(Converged, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, rc)
    
    if (Converged.AND.Nc<4) then
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
