module parquet_util
  
  use mpi_mod
  use global_parameter
  use math_mod

  implicit none

  character(len=30), parameter :: FaddB   = 'F+B'
  character(len=30), parameter :: FaddF   = 'F+F'
  character(len=30), parameter :: FsubF   = 'F-F'
  character(len=30), parameter :: MinusF  = '-F'
  character(len=30), parameter :: MinusB  = '-B'
  character(len=30), parameter :: Bosonic   = 'Bosonic'
  character(len=30), parameter :: Fermionic = 'Fermionic'

  ! ... variables specific to the parquet approach ...
  logical :: DGA=.false.  ! determine if the local fully irreducible vertex is
                          ! used. In case of DGA=.False., parquet approximation
                          ! is used.
  integer :: Nx           ! linear dimension of a1. 
  integer :: Ny           ! linear dimension of a2. 
  integer :: Nf           ! num of Matsubara frequencies
  integer :: Nt           ! linear dimension of combined momentum and frequency
                          ! variable 
  integer :: Nc           ! total number of sites
  integer :: Nb           ! the number of bosonic variables used for the massive
                          ! parallization.
  real(dp) :: f_damping   ! damping factor on the reducible vertex

  complex,     allocatable :: Sigma_H(:, :)  ! Hatree energy

  real(dp),    allocatable :: Ek(:, :)       ! tight-binding dispersion
  
  complex(dp), allocatable :: Delta(:)       ! DMFT hybridization function
  
  complex(dp), allocatable :: L_d(:, :, :)   ! local fully irreducible vertex in
                                             ! the density channel
  complex(dp), allocatable :: L_m(:, :, :)   !  ... in the magentic channel
  complex(dp), allocatable :: L_s(:, :, :)   !  ... in the singlet channel
  complex(dp), allocatable :: L_t(:, :, :)   !  ... in the triplet channel
  
  ! complete vertex in each channel
  complex(dp), allocatable :: F_d(:, :, :)     ! density channel
  complex(dp), allocatable :: F_m(:, :, :)     ! magnetic channel
  complex(dp), allocatable :: F_s(:, :, :)     ! singlet channel
  complex(dp), allocatable :: F_t(:, :, :)     ! triplet channel

  ! irreducible vertex in each channel 
  complex(dp), allocatable :: G_d(:, :, :)     ! density channel
  complex(dp), allocatable :: G_m(:, :, :)     ! magnetic channel
  complex(dp), allocatable :: G_s(:, :, :)     ! singlet channel
  complex(dp), allocatable :: G_t(:, :, :)     ! triplet channel

  ! single-particle property
  complex(dp), allocatable :: Gkw(:)           ! Green's function
  complex(dp), allocatable :: Sigma(:)         ! Self-energy

  complex(dp), allocatable :: chi0_ph(:)       ! bubble diagram in p-h channel
  complex(dp), allocatable :: chi0_pp(:)       ! bubble diagram in p-p channel

  complex(dp), allocatable :: k2_d1(:,:,:)     ! kernel approximation of Phi_d
  complex(dp), allocatable :: k2_m1(:,:,:)     ! kernel approximation of Phi_m
  complex(dp), allocatable :: k2_s1(:,:,:)     ! kernel approximation of Phi_s
  complex(dp), allocatable :: k2_t1(:,:,:)     ! kernel approximation of Phi_t
  complex(dp), allocatable :: k2_d2(:,:,:)     ! kernel approximation of Phi_d
  complex(dp), allocatable :: k2_m2(:,:,:)     ! kernel approximation of Phi_m
  complex(dp), allocatable :: k2_s2(:,:,:)     ! kernel approximation of Phi_s
  complex(dp), allocatable :: k2_t2(:,:,:)     ! kernel approximation of Phi_t

  complex(dp), allocatable :: mat(:,:,:)       ! temporary arry for parquet equation and self-energy calculations

  ! FFT 
  real(dp), allocatable :: C_wave_x(:)         ! work array for FFT in x
  real(dp), allocatable :: C_wave_y(:)         ! work array for FFT in y

  type :: indxmap
     integer :: ix
     integer :: iy
     integer :: iw
  end type indxmap

  type(indxmap), allocatable :: Index_fermionic(:), Index_bosonic(:)

contains
  !------------------------------------------------------------------------------
  subroutine readin
    implicit none
    
    !
    !  Purpose
    ! =========
    !   read in the parameters on each node. Default values are set in
    !   global_parameter.
    !   initialize the vertex functions, and prepare workspace for FFT.
    !
    !   There are three mapping conventions used in this program, they
    !   relate to each other as follows:
    !
    !   ix: the momentum/site index in x direction
    !   iy: the momentum/site index in y direction
    !   iw: the Matsubara frequency index
    !
    !   (ix, iy, iw) -> indx = [(ix-1)*Ny+iy-1]*Nf+iw
    !   
    !   ib: the index on each node for the assigned copies of the full vertex
    !
    !   indx = id*Nb + ib, where id is the node index between [0, ntasks)
    !   
    !   the fermionic and bosonic frequencies are given as
    !    
    !    v = Pi/beta * [2*(iw-Nf/2-1)+1] ;  v = Pi/beta*[-Nf+1 : Nf-1: 2]
    !    w = Pi/beta * [2*(iw-1)];          w = Pi/beta*[0: Nf-2 : 2]
    !

    ! ... local vars ...
    integer   :: i, j, k, itemp, idx
    real(dp)  :: t1, t2, t3, t4, dkx, dky, kx, ky
    
    open(unit=1, file='input.parquet', status='old')
    read(1, *) Nx, Ny, nOmega, beta, xU, nParticle, f_damping     
    read(1, *) DGA                              ! if the local fully irreducible
                                                ! vertex from DMFT is used or
                                                ! not.
    close(1)

    Nc = Nx*Ny
    Nf = 2*nOmega
    Nt = Nc*Nf
    Nb = Nt/ntasks/2
    
    if (mod(Nt, ntasks) /= 0) then
       call error_msg('The allocated total num. of CPU is not compatible for parallization.')
    else
       if (id == master) then
          write(*, '(a,i4,a)'), 'each node gets ', Nb, ' copies of Nc*Nf data.'
       end if
    end if

    if (.NOT. allocated(Delta)) allocate(Delta(Nf))
    if (Nc == 1) then
       open(unit=1, file='Delta.dat', status='unknown')
       do i = Nf/2+1, Nf
          read(1, *) t1, t2, t3
          Delta(i) = dcmplx(t2, t3)
          Delta(Nf-i+1) = conjg(Delta(i))
       end do
       close(1)
    else
       Delta = Zero
    end if

    !  --- set up the index map for the convenience of reference ---
    if (.NOT. allocated(Index_fermionic)) allocate(Index_fermionic(Nt))
    if (.NOT. allocated(Index_bosonic))   allocate(Index_bosonic(Nt/2))
    do i = 1, Nx
       do j = 1, Ny
          do k = 1, Nf
             idx = ((i-1)*Ny+j-1)*Nf + k
             Index_fermionic(idx) = indxmap(i, j, k) 
          end do
          do k = 1, Nf/2
             idx = ((i-1)*Ny+j-1)*Nf/2 + k
             Index_bosonic(idx) = indxmap(i, j, k)
          end do
       end do
    end do
 
    ! --- set up the tight-binding dispersion ---
    if (.NOT. allocated(Ek)) allocate(Ek(Nx, Ny))
    dkx = Two*Pi/Nx
    dky = Two*Pi/Ny
    Ek = Zero
    do i = 1, Nx
       kx = dkx*(i-1)
       do j = 1, Ny
          ky = dky*(j-1)
          if (Nx > 1) Ek(i, j) = Ek(i, j) - Two*cos(kx)
          if (Ny > 1) Ek(i, j) = Ek(i, j) - Two*cos(ky)
       end do
    end do
    if (Nc == 1) Ek = Zero

    if (.NOT. allocated(Sigma_H)) allocate(Sigma_H(Nx, Ny))
    Sigma_H = Zero
  
    if (.NOT. allocated(C_wave_x)) allocate(C_wave_x(4*Nx+15))
    if (.NOT. allocated(C_wave_y)) allocate(C_wave_y(4*Ny+15))   
    call Zffti(Nx, C_wave_x)
    call Zffti(Ny, C_wave_y)

    ! allocate arrays for the local fully irreducible vertex in each node
    if (.NOT. allocated(L_d)) allocate(L_d(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_m)) allocate(L_m(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_s)) allocate(L_s(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_t)) allocate(L_t(Nf, Nf, Nf/2))
    L_d =  xU
    L_m = -xU
    L_s = Two*xU
    L_t = Zero

    if (dga) then
       open(unit=1, file='L_ph.dat', status='old')
       do i = 1, Nf
          do j = 1, Nf
             do k = 1, Nf/2
                read(1, *) itemp, itemp, itemp, t1, t2, t3, t4
                L_d(i, j, k) = dcmplx(t1, t2)      ! density component  
                L_m(i, j, k) = dcmplx(t3, t4)      ! magnetic component 
             end do
          end do
       end do
       close(1)

       open(unit=1, file='L_pp.dat', status='old')
       do i = 1, Nf
          do j = 1, Nf
             do k = 1, Nf/2
                read(1, *) itemp, itemp, itemp, t1, t2, t3, t4
                L_s(i, j, k) = dcmplx(t1, t2)     ! singlet component 
                L_t(i, j, k) = dcmplx(t3, t4)     ! triplet component
             end do
          end do
       end do
       close(1)
    end if
    
    if (.NOT. allocated(G_d)) allocate(G_d(Nt, Nt, Nb))
    if (.NOT. allocated(G_m)) allocate(G_m(Nt, Nt, Nb))
    if (.NOT. allocated(G_s)) allocate(G_s(Nt, Nt, Nb))
    if (.NOT. allocated(G_t)) allocate(G_t(Nt, Nt, Nb))    
    G_d =  xU
    G_m = -xU
    G_s = Two*xU
    G_t = Zero

   ! initiate the complete vertex
    if (.NOT. allocated(F_d)) allocate(F_d(Nt, Nt, Nb))
    if (.NOT. allocated(F_m)) allocate(F_m(Nt, Nt, Nb))
    if (.NOT. allocated(F_s)) allocate(F_s(Nt, Nt, Nb))
    if (.NOT. allocated(F_t)) allocate(F_t(Nt, Nt, Nb))
    F_d = G_d
    F_m = G_m
    F_s = G_s
    F_t = G_t


    if (dga) then
       do k = 1, Nb
          ! corresponding component in the complete list is given at "id*Nb+k"
          if (mod(id*Nb+k, Nf/2) == 0 ) then
             idx = Nf/2
          else
             idx = mod(id*Nb+k, Nf/2)
          end if
          do i = 1, Nc
             do j = 1, Nc
                F_d(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = L_d(1:Nf, 1:Nf, idx)
                F_m(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = L_m(1:Nf, 1:Nf, idx)
                F_s(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = L_s(1:Nf, 1:Nf, idx)
                F_t(Nf*(i-1)+1:Nf*i, Nf*(j-1)+1:Nf*j, k) = L_t(1:Nf, 1:Nf, idx)
             end do
          end do
       end do
    end if

   end subroutine readin

  !------------------------------------------------------------------------------
  integer function list_index(P, typ) result(idx)
     
    type(Indxmap), intent(in) :: P
    character(len=30), intent(in) :: typ

    if (typ == 'Bosonic') then
       idx = ((P%ix-1)*Ny+P%iy-1)*Nf/2+P%iw
    elseif (typ == 'Fermionic') then
       idx =((P%ix-1)*Ny+P%iy-1)*Nf+P%iw
    else
       write(*, *) 'unknown Matsubara Frequency type!'
       stop
    end if

  end function list_index

  !------------------------------------------------------------------------------
  subroutine index_operation(idx1, idx2, operation, final_Indx)
    !
    ! Purpose
    ! =======
    !  For given two indices in the complete list, find out the resulting index for
    !  a given operation.
    !
    type(Indxmap), intent(in)     :: idx1, idx2
    character(len=30), intent(in) :: operation
    type(Indxmap), intent(out)    :: final_Indx

    ! ... local vars ...
    integer :: i, j, k

    if (Trim(operation) == FaddB) then 
       i = idx1%ix + idx2%ix - 1
       if (i > Nx) i = i - Nx
       j = idx1%iy + idx2%iy - 1
       if (j > Ny) j = j - Ny
       k = idx1%iw + idx2%iw - 1      
       final_Indx = indxmap(i, j, k)
    end if

    if (Trim(operation) == FaddF) then
       i = idx1%ix + idx2%ix - 1
       if (i > Nx) i = i - Nx
       j = idx1%iy + idx2%iy - 1
       if (j > Ny) j = j - Ny
       k = idx1%iw + idx2%iw - Nf 
       final_Indx = indxmap(i, j, k)
    end if

    if (Trim(operation) == MinusF) then
       i = -idx1%ix + Nx + 2
       if (i > Nx) i = i - Nx
       j = -idx1%iy + Ny + 2
       if (j > Ny) j = j - Ny
       k = -idx1%iw + Nf + 1
       final_Indx = indxmap(i, j, k)
    end if

    if (Trim(operation) == MinusB) then
       i = -idx1%ix + Nx + 2
       if (i > Nx) i = i - Nx
       j = -idx1%iy + Ny + 2
       if (j > Ny) j = j - Ny
       k = -idx1%iw + 2
       final_indx = indxmap(i, j, k)
    end if
  end subroutine index_operation

  !------------------------------------------------------------------------------
  subroutine pa_Gkw_Chi0(ite, Grt)
    !
    ! Purpose
    ! =======
    !   Generate the Green's function Gkw. The initial value is just the non-interacting
    !   one. The bubble diagrams of two convoluted Green's function are then 
    !   calculated from the Fourier-Transformed Gkw, i.e. Grt. The inverse Fourier
    !   transform is carried out by using the supplementation of spline interpolation.
    !
    ! Notation
    ! ========
    !   fermionic frequency: k=[1, Nf] -> v = Pi/beta* (Two*(k-Nf/2-1) + One)
    !   bosonic   frequency: k=[1, Nf] -> w = Pi/beta* Two*(k-Nf/2)
    !
    ! Returned values
    ! ===============
    !   Gkw and Chi0_ph, Chi0_pp will be updated
    !
    integer,     intent(in)  :: ite          ! loop index
    complex(dp), intent(out) :: Grt(Nx, Ny, Nf)

    ! ... local vars ...
    integer     :: i, j, k, idx, iTau
    integer     :: i1, j1
    real(dp)    :: t1, mu_UpperBound, mu_LowerBound
    real(dp)    :: dkx, dky, kx, ky, w, tau, dummy, FD1, FD2
    complex(dp) :: dummy1D(Nf), coutdata(Nf/2)
    complex(dp) :: Chi0rt_ph(Nx, Ny, Nf), Chi0rt_pp(Nx, Ny, Nf)
    character(len=30) :: FLE, str1
    character(len=10) :: Mtype


    if (.NOT. allocated(Gkw))     allocate(Gkw(Nt))
    if (.NOT. allocated(Sigma))   allocate(Sigma(Nt))
    if (ite == 1) Sigma = Zero
    do i = 1, Nx
       do j = 1, Ny
          idx = (Ny*(i-1)+j)*Nf
          Sigma_H(i, j) = Sigma(idx)
          if (id == master) write(*, "('Hartree energy for', 2i4, ' is', 2f12.6 )") i, j, Sigma_H(i, j)
       end do
    end do    

    ! --- first to adjust the chemical potential ---
    mu_UpperBound = mu + xU
    mu_LowerBound = mu - xU
    mu = Zero

111 continue
    dkx = Two*Pi/Nx
    dky = Two*Pi/Ny
    do i = 1, Nx
       do j = 1, Ny
          do k = 1, Nf
             w   = Pi/beta*(Two*(k-Nf/2-1)+One)
             idx = (Ny*(i-1)+j-1)*Nf + k
             if (Nc == 1) then
                Gkw(k) = One/(xi*w + mu - Delta(k) - Sigma(k))
             else                
                Gkw(idx) = One/(xi*w + mu - Ek(i, j) - Sigma(idx))
             end if
             dummy1D(k) = (Gkw(idx) - One/(xi*w + mu - Ek(i, j)))  ! Ek is zero when Nc = 1
          end do
          ! Fourier transform from w -> tau space with the supplemented
          ! function (non-interacting Green's function).
          do iTau = 1, Nf
             tau = beta/(Nf-1)*(iTau-1)
             dummy = Zero
             do k = 1, Nf
                w = Pi/beta*(Two*(k-Nf/2)-One)
                dummy = dummy + dble( exp(-xi*w*tau)*dummy1D(k) )/beta
             end do
             Grt(i, j, iTau) = dummy - exp((beta-tau)*(Ek(i, j)-mu))/(One + exp(beta*(Ek(i, j)-mu)))
          end do
       end do
    end do

    ! FFT from k to r space.
    do k = 1, Nf
       call fftf2d(Nx, Ny, Grt(1:Nx, 1:Ny, k), C_wave_x, C_wave_y)
    end do

    t1 = -Two*Grt(1, 1, Nf)

    if (id == master) write(*, "(' particle number is:', f12.6, ' chemical potential is', f12.6)") t1, mu
    if (abs(t1 - nParticle) > 1.d-4) then
       if (t1 > nParticle) then
          mu_UpperBound = mu
       else
          mu_LowerBound = mu
       end if
       mu = Half*(mu_UpperBound + mu_LowerBound)
       goto 111
    end if
   
    ! output the Green's function in k-w space
    if (id == master) then
       write(str1, '(I0.3)') ite
       FLE = 'Gkw-'//trim(str1)//'.dat'
       open(unit=1, file=FLE, status='unknown')    
       do i = 1, Nx
          kx = dkx*(i-1)
          do j = 1, Ny
             ky = dky*(j-1)
             do k = 1, Nf
                idx = (Ny*(i-1)+j-1)*Nf + k
                write(1, '(2f12.6, 3f20.12)') kx, ky, Pi/beta*(Two*(k-Nf/2)-One), Gkw(idx)
             end do
             write(1, *)
          end do
       end do
       close(1)
       
       ! output the Green's function in r-t space
       FLE = 'Grt-'//trim(str1)//'.dat'
       open(unit=1, file=FLE, status='unknown')
       do i = 1, Nx
          do j = 1, Ny
             do k = 1, Nf
                write(1, '(2i4, 3f20.12)') i, j, beta/(Nf-1)*(k-1), Grt(i, j, k)
             end do
             write(1, *)
          end do
       end do
       close(1)
    end if

    !
    ! Now, determine the bubble diagram of two convoluted Green's function, which will be 
    ! used to approximate the reducible vertex in each channel for components that are not 
    ! avaliable from the updated complete- and irreducible-vertex in this channel.
    !

    if (.NOT. allocated(Chi0_ph)) allocate(Chi0_ph(Nt/2))
    if (.NOT. allocated(Chi0_pp)) allocate(Chi0_pp(Nt/2))

    ! --- particle-hole bubble ---
    do k = 1, Nf
       do i = 1, Nx
          if (i == 1) then
             i1 = i
          else
             i1 = Nx-i+2
          end if
          do j = 1, Ny
             if (j == 1) then
                j1 = 1
             else
                j1 = Ny-j+2
             end if
             Chi0rt_ph(i, j, k) = -Grt(i, j, k)*Grt(i1, j1, Nf-k+1)     
          end do
       end do
       call fftb2d(Nx, Ny, Chi0rt_ph(1:Nx, 1:Ny, k), C_wave_x, C_wave_y)
    end do

    Mtype = 'Bosonic'
    do i = 1, Nx
       do j = 1, Ny
          do k = 1, Nf
             dummy1D(k) = Chi0rt_ph(i, j, k) 
          end do
          call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
          call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
          do k = 1, Nf/2
             idx = (Ny*(i-1)+j-1)*Nf/2 + k
             Chi0_ph(idx) = coutdata(k)
          end do
       end do
    end do
    
    ! --- particle-particle bubble ---
    do k = 1, Nf
       do i = 1, Nx
          do j = 1, Ny
             Chi0rt_pp(i, j, k) = -Half*Grt(i, j, k)**2
          end do
       end do
       call fftb2d(Nx, Ny, Chi0rt_pp(1:Nx, 1:Ny, k), C_wave_x, C_wave_y)
    end do

    do i = 1, Nx
       do j = 1, Ny
          do k = 1, Nf
             dummy1D(k) = Chi0rt_pp(i, j, k)    
          end do
          call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
          call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
          do k = 1, Nf/2
             idx = (Ny*(i-1)+j-1)*Nf/2 + k
             Chi0_pp(idx) = coutdata(k)
          end do          
       end do
    end do
    
    ! --- monitoring the bubble results ---
    if (id == master) then
       FLE = 'Chi0-'//trim(str1)//'.dat'
       open(unit=1, file=FLE, status='unknown')
       do i = 1, Nx
          kx = dkx*(i-1)
          do j = 1, Ny
             ky = dky*(j-1)
             do k = 1, Nf/2
                idx = (Ny*(i-1)+j-1)*Nf/2 + k
                write(1, '(2f12.6, 5f12.6)') kx, ky, Pi/beta*Two*(k-1), &
                     Chi0_ph(idx), Chi0_pp(idx) 
             end do
             write(1, *)
          end do
       end do
       close(1)
    end if

  end subroutine pa_Gkw_Chi0

  !------------------------------------------------------------------------------
  subroutine Memory_Release

    if (allocated(Delta)) Deallocate(Delta)

    ! clear the dispersion 
    if (allocated(Ek)) Deallocate(Ek)

    ! deallocate the Hatree energy
    if (allocated(Sigma_H)) Deallocate(Sigma_H)

    ! deallocate the index array
    if (allocated(Index_bosonic)) Deallocate(Index_bosonic)
    if (allocated(Index_fermionic)) Deallocate(Index_fermionic)

    ! deallocate arries for FFT
    if (allocated(C_wave_x)) Deallocate(C_wave_x)
    if (allocated(C_wave_y)) Deallocate(C_wave_y)    

    ! deallocate the local fully irreducible vertex 
    if (allocated(L_d)) Deallocate(L_d)
    if (allocated(L_m)) Deallocate(L_m)
    if (allocated(L_s)) Deallocate(L_s)
    if (allocated(L_t)) Deallocate(L_t)

    ! deallocate the complete veretx
    if (allocated(F_d)) Deallocate(F_d)
    if (allocated(F_m)) Deallocate(F_m)
    if (allocated(F_s)) Deallocate(F_s)
    if (allocated(F_t)) Deallocate(F_t) 

    ! deallocate the irreducible veretx in each channel
    if (allocated(G_d)) Deallocate(G_d)
    if (allocated(G_m)) Deallocate(G_m)
    if (allocated(G_s)) Deallocate(G_s)
    if (allocated(G_t)) Deallocate(G_t) 

    ! deallocate the single-particle Green's function
    if (allocated(Gkw))     Deallocate(Gkw)
    if (allocated(Sigma))   Deallocate(Sigma)
    if (allocated(Chi0_ph)) Deallocate(Chi0_ph)
    if (allocated(Chi0_pp)) Deallocate(Chi0_pp)


    if (allocated(K2_d1))  Deallocate(K2_d1)
    if (allocated(K2_m1))  Deallocate(k2_m1)
    if (allocated(K2_s1))  Deallocate(k2_s1)
    if (allocated(K2_t1))  Deallocate(k2_t1)
    if (allocated(K2_d2))  Deallocate(K2_d2)
    if (allocated(K2_m2))  Deallocate(k2_m2)
    if (allocated(K2_s2))  Deallocate(k2_s2)
    if (allocated(K2_t2))  Deallocate(k2_t2)

    ! deallocate temparary array
    if (allocated(mat))   Deallocate(mat)

  end subroutine Memory_Release

end module parquet_util
