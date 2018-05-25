module Parquet_kernel
  !
  ! Purpose
  ! ========
  !   determine the reducible vertex function in particle-hole and particle-particle channels
  ! 
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util

contains
  !-------------------------------------------------------------------------------------------
  subroutine reducible_vertex(ite)
     
    integer, intent(in) :: ite

    ! ... local vars ...
    integer     :: i, j, k, idx, idx1
    complex(dp) :: dummy, G1
    character(len=30) :: FLE, str
    
    complex(dp), allocatable :: dummy3D_1(:, :), dummy3D_2(:, :), dummy3D_3(:, :), dummy3D_4(:, :)
    type(Indxmap) :: ComIdx1, ComIdx2
    
    if (.NOT. allocated(dummy3D_1)) allocate(dummy3D_1(Nt, Nt))
    if (.NOT. allocated(dummy3D_2)) allocate(dummy3D_2(Nt, Nt))
    if (.NOT. allocated(dummy3D_3)) allocate(dummy3D_3(Nt, Nt))
    if (.NOT. allocated(dummy3D_4)) allocate(dummy3D_4(Nt, Nt))
    
    do k = 1, Nb
       idx = id*Nb + k
      
       ! start to calculate the reducible vertex function in particle-hole channel
       G1        = Zero_c
       dummy3D_3 = Zero_c
       dummy3D_4 = Zero_c

       ! --- density (d) and magnetic (m) channels ---
       !  Phi = Gamma *G*G* F
       do i = 1, Nt
          call index_operation(Index_fermionic(i), Index_bosonic(idx), FaddB, ComIdx1) ! k+q
          idx1 = list_index(ComIdx1, Fermionic)
          
          if (ComIdx1%iw > Nf .or. ComIdx1%iw < 1) then
             ! use the non-interacting Green's function when k+q is outside of the box
             dummy = One/(xi*Pi/beta*(Two*(ComIdx1%iw-Nf/2-1) + One) + mu - Ek(ComIdx1%ix, ComIdx1%iy))
          else
             dummy = Gkw(idx1) 
          end if
          do j = 1, Nt
             dummy3D_1(i, j) = F_d(i, j, k)*Gkw(i)*dummy/beta/Nc
             dummy3D_2(i, j) = F_m(i, j, k)*Gkw(i)*dummy/beta/Nc 
          end do
          G1 = G1 + xU*Gkw(i)*dummy*xU/beta/Nc
       end do
       
       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, G_d(1:Nt, 1:Nt, k), Nt, dummy3D_1(1:Nt, 1:Nt), Nt, Zero_c, dummy3D_3, Nt)
       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, G_m(1:Nt, 1:Nt, k), Nt, dummy3D_2(1:Nt, 1:Nt), Nt, Zero_c, dummy3D_4, Nt)

       ! Phi = F *G*G* Gamma
       do i = 1, Nt
          call index_operation(Index_fermionic(i), Index_bosonic(idx), FaddB, ComIdx1) ! k+q
          idx1 = list_index(ComIdx1, Fermionic)

          if (ComIdx1%iw > Nf .or. ComIdx1%iw < 1) then
             ! use the non-interacting Green's function when k+q is outside of the box
             dummy = One/(xi*Pi/beta*(Two*(ComIdx1%iw-Nf/2-1) + One) + mu - Ek(ComIdx1%ix, ComIdx1%iy))
          else
             dummy = Gkw(idx1)
          end if
          do j = 1, Nt
             dummy3D_1(j, i) = F_d(j, i, k)*Gkw(i)*dummy/beta/Nc
             dummy3D_2(j, i) = F_m(j, i, k)*Gkw(i)*dummy/beta/Nc
          end do
       end do

       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, dummy3D_1(1:Nt, 1:Nt), Nt, G_d(1:Nt, 1:Nt, k), Nt, One_c, dummy3D_3, Nt)
       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, dummy3D_2(1:Nt, 1:Nt), Nt, G_m(1:Nt, 1:Nt, k), Nt, One_c, dummy3D_4, Nt)
       
       do i = 1, Nt
          do j = 1, Nt
             if (ite == 1) then
                G_d(i, j, k) = (One-f_damping)*(xU*Chi0_ph(idx)*xU) + f_damping*(dummy3D_3(i, j) - G1 + xU*Chi0_ph(idx)*xU)
                G_m(i, j, k) = (One-f_damping)*(xU*Chi0_ph(idx)*xU) + f_damping*(dummy3D_4(i, j) - G1 + xU*Chi0_ph(idx)*xU)
             else
                G_d(i, j, k) = (One-f_damping)*(F_d(i, j, k)-G_d(i, j, k)) + f_damping*(dummy3D_3(i, j) - G1 + xU*Chi0_ph(idx)*xU)
                G_m(i, j, k) = (One-f_damping)*(F_m(i, j, k)-G_m(i, j, k)) + f_damping*(dummy3D_4(i, j) - G1 + xU*Chi0_ph(idx)*xU)
             end if
          end do
       end do
    
       ! --- singlet (s) and triplet (t) channel --- 
       G1        = Zero_c
       dummy3D_3 = Zero_c
       dummy3D_4 = Zero_c

       ! Psi = Gamma *G*G* F
       do i = 1, Nt
          call index_operation(Index_fermionic(i), Index_fermionic(i), MinusF, ComIdx1)    !  -k
          call index_operation(ComIdx1, Index_bosonic(idx), FaddB, ComIdx2)                ! q-k
          idx1 = list_index(ComIdx2, Fermionic)
          if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
             ! use the non-interacting green's function when q-k is outside the box
             dummy = One/( xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) )
          else
             dummy = Gkw(idx1)
          end if
          do j = 1, Nt
             dummy3D_1(i, j) = -Half*F_s(i, j, k)*Gkw(i)*dummy/beta/Nc
             dummy3D_2(i, j) =  Half*F_t(i, j, k)*Gkw(i)*dummy/beta/Nc
          end do
          G1 = G1 - xU*Gkw(i)*dummy*(Two*xU)/beta/Nc
       end do

       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, G_s(1:Nt, 1:Nt, k), Nt, dummy3D_1(1:Nt, 1:Nt), Nt, Zero_c, dummy3D_3, Nt)
       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, G_t(1:Nt, 1:Nt, k), Nt, dummy3D_2(1:Nt, 1:Nt), Nt, Zero_c, dummy3D_4, Nt)

       ! Psi = F *G*G* Gamma
       do i = 1, Nt
          call index_operation(Index_fermionic(i), Index_fermionic(i), MinusF, ComIdx1)    !  -k
          call index_operation(ComIdx1, Index_bosonic(idx), FaddB, ComIdx2)                ! q-k
          idx1 = list_index(ComIdx2, Fermionic)
          if (ComIdx2%iw > Nf .or. ComIdx2%iw < 1) then
             ! use the non-interacting green's function when q-k is outside the box
             dummy = One/( xi*Pi/beta*(Two*(ComIdx2%iw-Nf/2-1) + One) + mu - Ek(ComIdx2%ix, ComIdx2%iy) )
          else
             dummy = Gkw(idx1)
          end if
          do j = 1, Nt
             dummy3D_1(j, i) = -Half*F_s(j, i, k)*Gkw(i)*dummy/beta/Nc
             dummy3D_2(j, i) =  Half*F_t(j, i, k)*Gkw(i)*dummy/beta/Nc
          end do
       end do

       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, dummy3D_1(1:Nt, 1:Nt), Nt, G_s(1:Nt, 1:Nt, k), Nt, One_c, dummy3D_3, Nt)
       call ZGEMM('N', 'N', Nt, Nt, Nt, Half_c, dummy3D_2(1:Nt, 1:Nt), Nt, G_t(1:Nt, 1:Nt, k), Nt, One_c, dummy3D_4, Nt)

       do i = 1, Nt
          do j = 1, Nt
             if (ite == 1) then
                G_s(i, j, k) = (One-f_damping)*(Two*xU*Chi0_pp(idx)*Two*xU) + f_damping*(dummy3D_3(i, j) - G1 + (Two*xU)*Chi0_pp(idx)*(Two*xU))
                G_t(i, j, k) = f_damping*dummy3D_4(i, j)
             else
                G_s(i, j, k) = (One-f_damping)*(F_s(i, j, k)-G_s(i, j, k)) + f_damping*(dummy3D_3(i, j) - G1 + (Two*xU)*Chi0_pp(idx)*(Two*xU))
                G_t(i, j, k) = (One-f_damping)*(F_t(i, j, k)-G_t(i, j, k)) + f_damping*dummy3D_4(i, j)
             end if
          end do
       end do
       
    end do

    write(str, '(I0.3,a,I0.3)') ite, '-', id
    if (Nc <= 2) then
       FLE = 'Reducible_V-'//trim(str)
       open(unit=1, file=FLE, status='unknown')
       do i = 1, Nt
          do j = 1, Nt
             do k = 1, Nb
                write(1, '(3i5, 8f20.12)') i, j, k, G_d(i, j, k), G_m(i, j, k), G_s(i, j, k), G_t(i, j, k)
             end do
          end do
       end do
    end if
    close(1)

    if (allocated(dummy3D_1)) deallocate(dummy3D_1)
    if (allocated(dummy3D_2)) deallocate(dummy3D_2)
    if (allocated(dummy3D_3)) deallocate(dummy3D_3)
    if (allocated(dummy3D_4)) deallocate(dummy3D_4)

    !!! G_d, G_m, G_s, G_t are now updated, which temporarily store the reducible vertex functions.
       
  end subroutine reducible_vertex

  ! ----------------------------------------------------------------------------------------------------------------
  subroutine get_kernel_function(ite)
    !
    ! Purpose
    ! =======
    !   From the exact reducible vertex function calculated in routine "reducible_Vertex" determine  
    ! the kernel functions by scanning the edge of it. 
    !
    integer, intent(in) :: ite

    ! ... local vars ...
    integer     :: ichannel, i, ic, jc, k
    complex(dp), allocatable :: dummy1(:,:,:), dummy2(:, :, :)
    character(len=30) :: FLE, str


    if (.NOT. allocated(K2_d1)) allocate(K2_d1(Nt, Nc, Nt/2))
    if (.NOT. allocated(K2_m1)) allocate(K2_m1(Nt, Nc, Nt/2))
    if (.NOT. allocated(k2_s1)) allocate(K2_s1(Nt, Nc, Nt/2))
    if (.NOT. allocated(K2_t1)) allocate(K2_t1(Nt, Nc, Nt/2))
    if (.NOT. allocated(K2_d2)) allocate(K2_d2(Nc, Nt, Nt/2))
    if (.NOT. allocated(K2_m2)) allocate(K2_m2(Nc, Nt, Nt/2))
    if (.NOT. allocated(K2_s2)) allocate(K2_s2(Nc, Nt, Nt/2))
    if (.NOT. allocated(K2_t2)) allocate(K2_t2(Nc, Nt, Nt/2))

    ! assign trivial initial values to kernel functions 
    if (ite == 1) then
       K2_d1 = Zero_c
       k2_d2 = Zero_c
       K2_m1 = Zero_c
       K2_m2 = Zero_c
       K2_s1 = Zero_c
       K2_s2 = Zero_c
       K2_t1 = Zero_c
       K2_t2 = Zero_c
    end if



    ! determine K2 on the master node and distribute it to other nodes
    if (.NOT. allocated(dummy1)) allocate(dummy1(Nt, Nc, Nb))
    if (.NOT. allocated(dummy2)) allocate(dummy2(Nc, Nt, Nb))

    dummy1 = Zero_c
    dummy2 = Zero_c
    do ichannel = 1, 4
       do k = 1, Nb
             ! here we will simply scan the edge of the reducible vertex function to get kernel
             do i = 1, Nt
                do jc = 1, Nc
                   select case (ichannel)
                   case (1)
                      dummy1(i, jc, k) = G_d(i, jc*Nf, k)
                      dummy2(jc, i, k) = G_d(jc*Nf, i, k)
                   case (2)
                      dummy1(i, jc, k) = G_m(i, jc*Nf, k)
                      dummy2(jc, i, k) = G_m(jc*Nf, i, k)
                   case (3)
                      dummy1(i, jc, k) = G_s(i, (jc-1)*Nf+1, k)
                      dummy2(jc, i, k) = G_s((jc-1)*Nf+1, i, k)
                   case (4)
                      dummy1(i, jc, k) = G_t(i, (jc-1)*Nf+1, k)
                      dummy2(jc, i, k) = G_t((jc-1)*Nf+1, i, k)
                   end select
                end do
             end do
       end do

       if (id == master) then
          select case (ichannel)
          case (1)
             K2_d1(1:Nt, 1:Nc, 1:Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
             K2_d2(1:Nc, 1:Nt, 1:Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
          case (2) 
             K2_m1(1:Nt, 1:Nc, 1:Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
             K2_m2(1:Nc, 1:Nt, 1:Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
          case (3)
             K2_s1(1:Nt, 1:Nc, 1:Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
             K2_s2(1:Nc, 1:Nt, 1:Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
          case (4)
             K2_t1(1:Nt, 1:Nc, 1:Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
             K2_t2(1:Nc, 1:Nt, 1:Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
          end select
          do i = 1, ntasks-1
             call MPI_RECV(dummy1, Nt*Nc*Nb, MPI_DOUBLE_COMPLEX, i, i, MPI_COMM_WORLD, status, rc)
             call MPI_RECV(dummy2, Nt*Nc*Nb, MPI_DOUBLE_COMPLEX, i, i, MPI_COMM_WORLD, status, rc)
             select case (ichannel)
             case (1)
                K2_d1(1:Nt, 1:Nc, (i*Nb+1):(i+1)*Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
                K2_d2(1:Nc, 1:Nt, (i*Nb+1):(i+1)*Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
             case (2)
                K2_m1(1:Nt, 1:Nc, (i*Nb+1):(i+1)*Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
                K2_m2(1:Nc, 1:Nt, (i*Nb+1):(i+1)*Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
             case (3)
                K2_s1(1:Nt, 1:Nc, (i*Nb+1):(i+1)*Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
                K2_s2(1:Nc, 1:Nt, (i*Nb+1):(i+1)*Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
             case (4)
                K2_t1(1:Nt, 1:Nc, (i*Nb+1):(i+1)*Nb) = dummy1(1:Nt, 1:Nc, 1:Nb)
                K2_t2(1:Nc, 1:Nt, (i*Nb+1):(i+1)*Nb) = dummy2(1:Nc, 1:Nt, 1:Nb)
             end select
          end do
       else
          call MPI_SEND(dummy1, Nt*Nc*Nb, MPI_DOUBLE_COMPLEX, master, id, MPI_COMM_WORLD, rc)
          call MPI_SEND(dummy2, Nt*Nc*Nb, MPI_DOUBLE_COMPLEX, master, id, MPI_COMM_WORLD, rc)
       end if
    end do

    call MPI_BCAST(K2_d1, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_d2, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_m1, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_m2, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_s1, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_s2, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_t1, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_t2, Nc*Nt*Nt/2, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)

    if (id == master) then
       write(str, '(I0.3)') ite
       FLE = 'K2-'//trim(str)//'.dat'
       open(unit = 1, file = FLE, status = 'unknown')
       do i = 1, Nt
          do ic = 1, Nc 
             do k = 1, Nt/2
                write(1, '(3i4, 16f10.5)') i, ic, k, K2_d1(i,ic,k), K2_d2(ic,i,k), K2_m1(i,ic,k), K2_m2(ic,i,k), &
                     K2_s1(i,ic,k), K2_s2(ic,i,k), K2_t1(i,ic,k), K2_t2(ic,i,k)
             end do
             write(1, *) 
          end do
       end do
       close(1)
    end if
    
    if (allocated(dummy1)) deallocate(dummy1)
    if (allocated(dummy2)) deallocate(dummy2)

  end subroutine get_kernel_function
 
  !-------------------------------------------------------------------------------------------------
  complex(dp) function kernel(channel, k1, k2, q)
    
    character(len=1), intent(in) :: channel
    type(Indxmap),    intent(in) :: k1, k2, q
  
    ! ... local vars ...
    integer       :: ic, jc, qc, iw, jw, qw, idx, idx1 
    type(Indxmap) :: k1_, k2_, q_
    complex(dp)   :: background
    
    if (q%iw > 0) then
       k1_ = k1
       k2_ = k2
        q_ = q
    else
       call index_operation(k1, k1, MinusF, K1_)
       call index_operation(k2, k2, MinusF, K2_)
       call index_operation(q, q, MinusB, q_)
    end if

    ic = (k1_%ix-1)*Ny+k1_%iy
    jc = (k2_%ix-1)*Ny+k2_%iy
    qc = (q_%ix-1)*Ny+q_%iy
    iw = k1_%iw
    jw = k2_%iw
    qw = q_%iw

    idx  = ((k1_%ix-1)*Ny+k1_%iy)*Nf
    idx1 = ((k1_%ix-1)*Ny+k1_%iy-1)*Nf+1 

    kernel = Zero
    if (qw > Nf/2) return

    select case (channel)
    case ('d')
    
          background = K2_d1(idx, jc, list_index(q_, Bosonic))

       if (iw > Nf .or. iw < 1) then  ! k1 is out of range
          if (jw > Nf .or. jw < 1) then  ! k2 if out of range
             kernel = background
          else
             kernel = K2_d2(ic, list_index(K2_, Fermionic), list_index(q_, Bosonic))
          end if
       else
          if (jw > Nf .or. jw < 1) then
             kernel = K2_d1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic))
          else
             Kernel = K2_d1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic)) + K2_d2(ic, list_index(k2_, Fermionic), list_index(q_, Bosonic)) - background
          end if
       end if       
    case ('m')
          background = K2_m1(idx, jc, list_index(q_, Bosonic))
       if (iw > Nf .or. iw < 1) then  ! k1 is out of range
          if (jw > Nf .or. jw < 1) then  ! k2 if out of range
             kernel = background
          else
             kernel = K2_m2(ic, list_index(K2_, Fermionic), list_index(q_, Bosonic))
          end if
       else
          if (jw > Nf .or. jw < 1) then
             kernel = K2_m1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic))
          else
             Kernel = K2_m1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic)) + K2_m2(ic, list_index(k2_, Fermionic), list_index(q_, Bosonic)) - background
          end if
       end if
    case ('s')
          background = K2_s1(idx1, jc, list_index(q_, Bosonic))
       if (iw > Nf .or. iw < 1) then  ! k1 is out of range
          if (jw > Nf .or. jw < 1) then  ! k2 if out of range
             kernel = background
          else
             kernel = K2_s2(ic, list_index(K2_, Fermionic), list_index(q_, Bosonic))
          end if
       else
          if (jw > Nf .or. jw < 1) then
             kernel = K2_s1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic))
          else
             Kernel = K2_s1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic)) + K2_s2(ic, list_index(k2_, Fermionic), list_index(q_, Bosonic)) - background
          end if
       end if
    case ('t')
          background = K2_t1(idx1, jc, list_index(q_, Bosonic))
       if (iw > Nf .or. iw < 1) then  ! k1 is out of range
          if (jw > Nf .or. jw < 1) then  ! k2 if out of range
             kernel = background
          else
             kernel = K2_t2(ic, list_index(K2_, Fermionic), list_index(q_, Bosonic))
          end if
       else
          if (jw > Nf .or. jw < 1) then
             kernel = K2_t1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic))
          else
             Kernel = K2_t1(list_index(k1_, Fermionic), jc, list_index(q_, Bosonic)) + K2_t2(ic, list_index(k2_, Fermionic), list_index(q_, Bosonic)) - background
          end if
       end if
    end select
     
    if (q%iw < 1)  Kernel = conjg(Kernel)

  end function kernel 
  

  
end module Parquet_kernel
