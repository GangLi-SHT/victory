module parquet_equation
  !
  ! Purpose
  ! ========
  !   in this module, there is only one routine which solves the parquet equations by using 
  ! the kernel function derived in module 'pa_kernel'.
  ! 
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util
  use parquet_kernel
  
  contains
  ! ------------------------------------------------------------------------------------------
  subroutine solve_parquet_equation(ite)
    implicit none
    
    integer, intent(in) :: ite
    
    ! ... local vars ...
    !integer(KIND=MPI_ADDRESS_KIND) :: iadummy !!! NEEDED ONLY IF non-blocking MPI_ISEND is used !!!
    integer       :: idx
    integer       :: i, j, k, i1, j1, k1, is
    integer       :: inode, ichannel
    type(Indxmap) :: ComIdx1, ComIdx2, ComIdx3, ComIdx4, ComIdx5, ComIdx2_s, ComIdx3_s
    type(Indxmap) :: map_i, map_k1, map_si, map_sk1! combined index for momentum and frequency
    character(len=10) :: FLE, str
    logical :: symm_list(Ns)
    
    if (.NOT. allocated(mat)) allocate(mat(Nt, Nt, Nb))
    
    ! --- kernel vertex + fully irreducible vertex
    F_d = G_d
    F_m = G_m
    F_s = G_s
    F_t = G_t
    do k = 1, Nb
      idx = index_bosonic_IBZ(id*Nb+k)%iw
      do i = 1, Nc
        do j = 1, Nc
          F_d(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = F_d(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) + L_d(1:Nf, 1:Nf, idx)
          F_m(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = F_m(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) + L_m(1:Nf, 1:Nf, idx)
          F_s(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = F_s(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) + L_s(1:Nf, 1:Nf, idx)
          F_t(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = F_t(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) + L_t(1:Nf, 1:Nf, idx)
        end do
      end do
    end do
    
    ! Important note:
    ! ===============  
    ! The communication between nodes needed here is provided by BCAST
    ! Depending on the MPI version, the message size in BCAST is limited to ~2GB (2^31 bytes); 
    ! In case you need to pass larger messages, it is necessary to split them in smaller chunks.
    ! Alternatively, if your Nb >1 you may use more cores for paralelization.
    ! Newest versions of MPI may support larger messages; 
    ! in this case, please remove the following line of code.
    
    if (Nt*Nt*Nb> 134217728)  call error_msg('The MPI message size exceeds the 2^31 byte limit.')
    
    do ichannel = 1, 4       
      do inode = 0, ntasks-1
        
        select case (ichannel)         
          case(1)                       ! (1) density channel                                             
             mat = G_d                                                                                    
          case(2)                       ! (2) magnetic channel  
             mat = G_m                                                                                    !!!
          case(3)                       ! (3) singlet channel   
             mat = G_s                                                                                    !!!
          case(4)                       ! (4) triplet channel  
             mat = G_t                                                                                    
        end select                                                                                      
        
        call MPI_BCAST(mat, Nt*Nt*Nb, MPI_DOUBLE_COMPLEX, inode, MPI_COMM_WORLD, rc)                    
        
        select case (ichannel)
          case (1, 2) ! density and magnetic channel
          !  
          !  rotation 1: Phi(k, k+q; k'-k) -> Phi(k1, k1'; q1) 
          !
            do i = 1, Nt      ! k
              
              map_i = index_fermionic(i)
              call index_minusF(map_i, ComIdx5) ! -k            
              
              do k = 1, Nb  ! q                   
                do k1 = 1, Nb   ! q1
                  
                  map_k1= index_bosonic_IBZ(inode*Nb+k1)
                  call list_symmetries(lattice_type,map_k1,symm_list)
                  
                  !symmetries    
                  do is= 1, Ns 
                    if (symm_list(is)) then
                    ! is=1  
                      call index_minusB(map_k1, ComIdx1) ! -q1                 
                      call symmetry_operation(lattice_type,is,ComIdx1, ComIdx4) ! -S(q1)=S(-q1)
                      call symmetry_operation(lattice_type,is,map_k1, map_sk1) ! S(q1)
                      
                      call index_FaddB(map_i, index_bosonic_IBZ(id*Nb+k), ComIdx2)  ! k+q                    
                      call symmetry_operation_inv(lattice_type,is,map_i,ComIdx2,map_si,ComIdx2_s ) ! S^-1(k), S^-1(k+q)
                      
                      i1=list_index_F(map_si)                   
                      
                      call index_FaddB(map_i, map_sk1, ComIdx3) ! k'=k+S(q1) 
                      
                      j = list_index_F(ComIdx3)
                      j1 = list_index_F(ComIdx2_s) 
                      if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then !check if frequency index of k' is within box
                        if (ComIdx2_s%iw >=1 .and. ComIdx2_s%iw <= Nf) then !check if frequency index of S^-1(k+q) is within box
                          if (ichannel == 1) then
                            F_d(i, j, k) = F_d(i, j, k) - 0.5d0*mat(i1, j1, k1)
                            F_m(i, j, k) = F_m(i, j, k) - 0.5d0*mat(i1, j1, k1)
                          elseif (ichannel == 2) then
                            F_d(i, j, k) = F_d(i, j, k) - 1.5d0*mat(i1, j1, k1)
                            F_m(i, j, k) = F_m(i, j, k) + 0.5d0*mat(i1, j1, k1)
                          end if
                        else
                          if (ichannel == 1) then
                            F_d(i, j, k) = F_d(i, j, k) - 0.5d0*Kernel('d', map_si, ComIdx2_s, map_k1)
                            F_m(i, j, k) = F_m(i, j, k) - 0.5d0*Kernel('d', map_si, ComIdx2_s, map_k1)
                          elseif (ichannel == 2) then
                            F_d(i, j, k) = F_d(i, j, k) - 1.5d0*Kernel('m', map_si, ComIdx2_s, map_k1)
                            F_m(i, j, k) = F_m(i, j, k) + 0.5d0*Kernel('m', map_si, ComIdx2_s, map_k1)
                          end if
                        end if
                      end if
                      
                      if (map_k1%iw > 1) then !if omega non-zero              
                        
                        ! Time reversal symmetric part: Phi(k, k+q; k'-k) -> Phi(-k1, -k1'; -q1)=conjg(Phi(k1, k1'; q1))
                        call index_MinusF(ComIdx2, ComIdx3) ! -k-q
                        call symmetry_operation_inv(lattice_type,is,ComIdx5,ComIdx3,map_si,ComIdx2_s ) ! S^-1(-k), S^-1(-k-q)
                        
                        i1=list_index_F(map_si)         
                        
                        call index_FaddB(map_i, ComIdx4, ComIdx3) ! k'=k+S(-q1) 
                        
                        if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then !check if frequency index of k' is within box
                          j = list_index_F(ComIdx3)            
                          j1 = list_index_F(ComIdx2_s)   
                          
                          if (ComIdx2_s%iw >=1 .and. ComIdx2_s%iw <= Nf) then !check if frequency index of S^-1(-k-q) is within box
                            if (ichannel == 1) then
                              F_d(i, j, k) = F_d(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                              F_m(i, j, k) = F_m(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                            elseif (ichannel == 2) then
                              F_d(i, j, k) = F_d(i, j, k) - 1.5d0*conjg(mat(i1, j1, k1))
                              F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                            end if
                          else
                            if (ichannel == 1) then
                              F_d(i, j, k) = F_d(i, j, k) - 0.5d0*conjg(Kernel('d', map_si, ComIdx2_s, map_k1))
                              F_m(i, j, k) = F_m(i, j, k) - 0.5d0*conjg(Kernel('d', map_si, ComIdx2_s, map_k1))
                            elseif (ichannel == 2) then
                              F_d(i, j, k) = F_d(i, j, k) - 1.5d0*conjg(Kernel('m', map_si, ComIdx2_s, map_k1))
                              F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(Kernel('m', map_si, ComIdx2_s, map_k1))
                            end if
                          end if
                        end if               
                      end if !q1%iw>0
                    
                    end if !symm_list
                  end do !symmetries
                end do !Nb
              end do! Nb
            end do !Nt   
            !
            ! rotation 2: Phi(k, q-k'; k'-k) -> Phi(k1, k1'; q1)
            !
            do i = 1, Nt   ! k
              map_i = index_fermionic(i)
              call index_minusF(map_i, ComIdx5) ! -k            
              do k = 1, Nb  ! q                   
                do k1 = 1, Nb   ! q1
                  
                  map_k1= index_bosonic_IBZ(inode*Nb+k1)
                  call list_symmetries(lattice_type,map_k1,symm_list)
                  
                  !symmetries    
                  do is= 1, Ns 
                    if (symm_list(is)) then
                      !  is=1
                      call index_minusB(map_k1, ComIdx1) ! -q1                 
                      call symmetry_operation(lattice_type,is,ComIdx1, ComIdx4) ! -S(q1)=S(-q1)
                      call symmetry_operation(lattice_type,is,map_k1, map_sk1) ! S(q1)
                      
                      call index_FaddB(ComIdx5, index_bosonic_IBZ(id*Nb+k), ComIdx2)  ! -k+q                    
                      call symmetry_operation_inv(lattice_type,is,map_i,ComIdx2,map_si,ComIdx2_s ) ! S^-1(k), S^-1(-k+q)
                      
                      i1=list_index_F(map_si)                   
                      
                      call index_FaddB(map_i, map_sk1, ComIdx3) ! k'=k+S(q1) 
                      call index_FaddB(ComIdx2_s, ComIdx1,  ComIdx3_s) ! k'=S^-1(q-k)-q1 
                      
                      j = list_index_F(ComIdx3)
                      j1 = list_index_F(ComIdx3_s)
                      
                      if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then !check if frequency index of k+S(q1) is within box
                        if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then !check if frequency index of S^-1(q-k)-q1 is within box
                          if (ichannel == 1) then
                            F_s(i, j, k) = F_s(i, j, k) + 0.5d0*mat(i1, j1, k1)
                            F_t(i, j, k) = F_t(i, j, k) - 0.5d0*mat(i1, j1, k1)
                          elseif (ichannel == 2) then
                            F_s(i, j, k) = F_s(i, j, k) - 1.5d0*mat(i1, j1, k1)
                            F_t(i, j, k) = F_t(i, j, k) - 0.5d0*mat(i1, j1, k1)
                          end if
                        else
                          if (ichannel == 1) then
                            F_s(i, j, k) = F_s(i, j, k) + 0.5d0*Kernel('d', map_si, ComIdx3_s, map_k1)
                            F_t(i, j, k) = F_t(i, j, k) - 0.5d0*Kernel('d', map_si, ComIdx3_s, map_k1)
                          else
                            F_s(i, j, k) = F_s(i, j, k) - 1.5d0*Kernel('m', map_si, ComIdx3_s, map_k1)
                            F_t(i, j, k) = F_t(i, j, k) - 0.5d0*Kernel('m', map_si, ComIdx3_s, map_k1)
                          end if
                        end if
                      end if
                      
                      if (map_k1%iw > 1) then !if omega non-zero
                        ! Time reversal symmetric part: Phi(k, q-k'; k'-k) -> Phi(-k1, -k1'; -q1)=conjg(Phi(k1, k1'; q1))
                        
                        call index_minusF(ComIdx2, ComIdx3) ! k-q
                        call symmetry_operation_inv(lattice_type,is,ComIdx5,ComIdx3,map_si,ComIdx2_s ) ! S^-1(-k), S^-1(k-q)
                        i1=list_index_F(map_si)         
                        call index_FaddB(map_i, ComIdx4, ComIdx3) ! k'=k+S(-q1) 
                        call index_FaddB(ComIdx2_s, ComIdx1, ComIdx3_s) ! -k'=S^-1(k-q)-q1                     
                        
                        if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then !check if frequency index of k+S(-q1) is within box
                          j = list_index_F(ComIdx3) 
                          j1 = list_index_F(ComIdx3_s)
                          if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then !check if frequency index of S^-1(k-q)-q1 is within box
                            if (ichannel == 1) then
                              F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                              F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                            elseif (ichannel == 2) then
                              F_s(i, j, k) = F_s(i, j, k) - 1.5d0*conjg(mat(i1, j1, k1))
                              F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                            end if
                          else
                            if (ichannel == 1) then
                              F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(Kernel('d', map_si, ComIdx3_s, map_k1))
                              F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(Kernel('d', map_si, ComIdx3_s, map_k1))
                            else
                              F_s(i, j, k) = F_s(i, j, k) - 1.5d0*conjg(Kernel('m', map_si, ComIdx3_s, map_k1))
                              F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(Kernel('m', map_si, ComIdx3_s, map_k1))
                            end if
                          end if                
                        end if
                      end if !q1%iw>0
                    end if !symm_list
                  end do !symmetries
                end do !Nb
              end do! Nb
            end do !Nt   
            !
            ! rotation 3: Phi(k, k'; q-k-k') -> Phi(k1, k1'; q1)
            !
            do i = 1, Nt   ! k
              
              map_i = index_fermionic(i)
              call index_minusF(map_i, ComIdx5) ! -k            
              
              do k = 1, Nb  ! q        
                do k1 = 1, Nb   ! q1
                  map_k1= index_bosonic_IBZ(inode*Nb+k1)
                  call list_symmetries(lattice_type,map_k1,symm_list)
                  
                  !symmetries    
                  do is= 1, Ns 
                    if (symm_list(is)) then
                      ! is=1
                      call index_minusB(map_k1, ComIdx1) ! -q1                 
                      call symmetry_operation(lattice_type,is,ComIdx1, ComIdx4) ! -S(q1)=S(-q1)
                      call symmetry_operation(lattice_type,is,map_k1, map_sk1) ! S(q1)
                      
                      call index_FaddB(ComIdx5, index_bosonic_IBZ(id*Nb+k),ComIdx2)  ! -k+q                    
                      call symmetry_operation_inv(lattice_type,is,map_i,ComIdx2,map_si,ComIdx2_s ) ! S^-1(k), S^-1(-k+q)
                      
                      i1=list_index_F(map_si)                   
                      
                      call index_FaddB(ComIdx2, ComIdx4, ComIdx3) ! k'=q-k-S(q1) 
                      call index_FaddB(ComIdx2_s, ComIdx1, ComIdx3_s) ! k'=S^-1(q-k)-q1 
                      
                      if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then !check if frequency index of q-k-S(q1) is within box
                        j = list_index_F(ComIdx3)
                          if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then !check if frequency index of S^-1(q-k)-q1 is within box
                          j1 = list_index_F(ComIdx3_s)
                          if (ichannel == 1) then
                            F_s(i, j, k) = F_s(i, j, k) + 0.5d0*mat(i1, j1, k1)
                            F_t(i, j, k) = F_t(i, j, k) + 0.5d0*mat(i1, j1, k1)
                          elseif (ichannel == 2) then
                            F_s(i, j, k) = F_s(i, j, k) - 1.5d0*mat(i1, j1, k1)
                            F_t(i, j, k) = F_t(i, j, k) + 0.5d0*mat(i1, j1, k1)
                          end if
                        end if
                      end if
                      
                      if (map_k1%iw > 1) then !if omega non-zero
                        ! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)   
                        call index_minusF(ComIdx2, ComIdx3) ! k-q
                        call symmetry_operation_inv(lattice_type,is,ComIdx5,ComIdx3,map_si,ComIdx2_s ) ! S^-1(-k), S^-1(k-q)
                        i1=list_index_F(map_si)         
                        call index_FaddB(ComIdx2, map_sk1,ComIdx3) ! k'=q-k+S(q1) 
                        call index_FaddB(ComIdx2_s, ComIdx1, ComIdx3_s) ! -k'=S^-1(k-q)-q1                     
                        
                        if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                          j = list_index_F(ComIdx3)                   
                          if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then
                            j1 = list_index_F(ComIdx3_s)
                            if (ichannel == 1) then
                              F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                              F_t(i, j, k) = F_t(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                            elseif (ichannel == 2) then
                              F_s(i, j, k) = F_s(i, j, k) - 1.5d0*conjg(mat(i1, j1, k1))
                              F_t(i, j, k) = F_t(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                            end if
                          end if
                        end if
                        
                      end if !q1%iw>0
                    end if !symm_list
                  end do !symmetries
                end do !Nb
              end do !Nb
            end do!Nt
            
          case (3, 4) ! spin and triplet channel
            !
            ! rotation 4: Phi(k, k'; k+k'+q) -> Phi(k1, k1'; q1) 
            !
            do i = 1, Nt   ! k
              
              map_i = index_fermionic(i)
              call index_minusF(map_i,  ComIdx5) ! -k
              ! i1 = i      ! k1 = k
              
              do k = 1, Nb   ! q
                
                call index_FaddB(map_i, index_bosonic_IBZ(id*Nb+k),  ComIdx1)  ! k+q
                call index_minusF(ComIdx1,  ComIdx2) ! -k-q
                
                do k1 = 1, Nb  ! q1   
                  map_k1= index_bosonic_IBZ(inode*Nb+k1)
                  call list_symmetries(lattice_type,map_k1,symm_list)
                  
                  !symmetries    
                  do is= 1, Ns 
                    if (symm_list(is)) then
                      
                      call symmetry_operation(lattice_type,is,map_k1, map_sk1) ! S(q1)
                      call symmetry_operation_inv(lattice_type,is,map_i,ComIdx2,map_si,ComIdx2_s ) ! S^-1(k), S^-1(-k-q)
                      
                      i1=list_index_F(map_si)
                      
                      call index_FaddB(ComIdx2, map_sk1, ComIdx3) ! k'=S(q1)-k-q 
                      call index_FaddB(ComIdx2_s, map_k1, ComIdx3_s) ! k'=q1+S^-1(-k-q) 
                      
                      if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then !check if frequency index of S(q1)-k-q is within box
                        
                        j = list_index_F(ComIdx3)
                        
                        if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then   
                          j1 = list_index_F(ComIdx3_s)
                          if (ichannel == 3) then
                            F_d(i, j, k) = F_d(i, j, k) + 0.5d0*mat(i1, j1, k1)
                            F_m(i, j, k) = F_m(i, j, k) - 0.5d0*mat(i1, j1, k1)
                          elseif (ichannel == 4) then
                            F_d(i, j, k) = F_d(i, j, k) + 1.5d0*mat(i1, j1, k1)
                            F_m(i, j, k) = F_m(i, j, k) + 0.5d0*mat(i1, j1, k1)
                          end if
                        end if
                      end if
                      
                      if (map_k1%iw > 1) then                        
                        ! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)   
                        call symmetry_operation_inv(lattice_type,is,ComIdx5,ComIdx1,map_si,ComIdx2_s ) ! S^-1(-k), S^-1(k+q)
                        i1=list_index_F(map_si)
                        
                        call index_minusB(map_sk1, ComIdx4) ! -S(q1)
                        call index_FaddB(ComIdx2, ComIdx4,  ComIdx3) ! k'=-S(q1)-k-q                    
                        
                        if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                          j = list_index_F(ComIdx3)              
                          call index_FaddB(ComIdx2_s, map_k1, ComIdx3_s) ! -k'=q1+S^-1(k+q)    
                          
                          if (ComIdx3_s%iw >=1 .and. ComIdx3_s%iw <= Nf) then
                            j1 = list_index_F(ComIdx3_s)
                            if (ichannel == 3) then
                              F_d(i, j, k) = F_d(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                              F_m(i, j, k) = F_m(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                            elseif (ichannel == 4) then
                              F_d(i, j, k) = F_d(i, j, k) + 1.5d0*conjg(mat(i1, j1, k1))
                              F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                            end if
                            
                          end if
                        end if
                      end if !q1%iw>0
                      
                    end if !symm_list
                  end do ! symmetries
                end do !Nb
              end do !Nb   
            end do !Nt
          
        end select
        
      end do
    end do
    
    !call MPI_Barrier(MPI_COMM_WORLD, rc)  ! syncretize all ranks !!! --- NEVER USE MPI_Barrier in production code --- !!!
    
    ! output complete vertex function
    write(str, '(I0.3,a1,I0.3)') ite,'-',id
    if (Nc == 1 ) then
      FLE = 'F-'//trim(str)
      open(unit=1, file=FLE, status='unknown')
      do i = 1, Nt
        do j = 1, Nt
          do k = 1, Nb
            write(1, '(3i5, 8f20.12)') i, j, k, F_d(i, j, k), F_m(i, j, k), F_s(i, j, k), F_t(i, j, k)
          end do
        end do
      end do
    end if
  close(1)        
  
  end subroutine solve_parquet_equation
  
end module parquet_equation
