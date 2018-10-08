module  parquet_util
  
  use mpi_mod
  use global_parameter
  use math_mod
  
  implicit none

  character(len=30), parameter :: Bosonic   = 'Bosonic'
  character(len=30), parameter :: Fermionic = 'Fermionic'
  character(len=30), parameter :: Square = 'Square'
  
  character(len=30) :: lattice_type
  
  ! ... variables specific to the parquet approach ...
  logical :: DGA=.false.  ! determine if the local fully irreducible vertex is
                          ! used. In case of DGA=.False., parquet approximation
                          ! is used.
  !logical :: Edge_Scan =.false.    ! scheme for determing the kernel function, .true. for 
                          ! scanning edge; .false. for calculation with previous
                          ! kernel functions. 
   
  logical :: Sigma_ini=.false.  ! if true, then the code uses user provided Sigma_ini
                          ! otherwise Sigma=0 for the first iteration
  logical :: No_sigma_update=.false.  ! if true, then at every iteration Sigma is set back to 
                            ! some input local Sigma (Sigma_DMFT)
  logical :: Eigen=.false.  ! if .true. eigenvalues are calculated at every iteration
  
  integer :: Nx           ! linear dimension of a1. 
  integer :: Ny           ! linear dimension of a2. 
  integer :: Nx_IBZ       ! linear dimension of the IBZ in x direction
  integer :: Ngrain       ! number of additonal k points in the coarse-grained k-sums
  
  integer :: Nf           ! num of Matsubara frequencies
  integer :: Nt           ! linear dimension of combined momentum and frequency
                          ! variable 
  integer :: Nred         ! linear dimension of combined momentum and frequency
                          ! variable, but reduced to the IBZ                      
  integer :: Ns           ! number of point group symmetries for a given lattice                        
  integer :: Nc           ! total number of sites
  integer :: Nb           ! the number of bosonic variables used for the massive
                          ! parallization.
  real(dp) :: f_damping   ! damping factor on the reducible vertex
  integer :: f_range      ! factor multiplying the range of frequencies out of the box for selfenergy
  
  real(dp) :: f_sigma     ! convergence condition for self_energy (relative error)
  real(dp) :: relative_err
  real(dp) :: eign_m, eign_s, eign_t, eign_pp
  
  integer :: ite_max    ! maximal number of iterations
  
  complex,     allocatable :: Sigma_H(:, :)  ! Hatree energy
  
  real(dp),    allocatable :: Ek(:, :)       ! tight-binding dispersion
  real(dp),    allocatable :: Ek_grain(:, :, :, :) ! tight-binding dispersion for a finer k grid
  
  complex(dp), allocatable :: Delta(:)       ! DMFT hybridization function
  
  complex(dp), allocatable :: L_d(:, :, :)   ! local fully irreducible vertex in
                                             ! the density channel
  complex(dp), allocatable :: L_m(:, :, :)   !  ... in the magentic channel
  complex(dp), allocatable :: L_s(:, :, :)   !  ... in the singlet channel
  complex(dp), allocatable :: L_t(:, :, :)   !  ... in the triplet channel
  
  !complete vertex in each channel
  complex(dp), allocatable :: F_d(:, :, :)     ! density channel
  complex(dp), allocatable :: F_m(:, :, :)     ! magnetic channel
  complex(dp), allocatable :: F_s(:, :, :)     ! singlet channel
  complex(dp), allocatable :: F_t(:, :, :)     ! triplet channel
  
  !irreducible vertex in each channel 
  complex(dp), allocatable :: G_d(:, :, :)     ! density channel
  complex(dp), allocatable :: G_m(:, :, :)     ! magnetic channel
  complex(dp), allocatable :: G_s(:, :, :)     ! singlet channel
  complex(dp), allocatable :: G_t(:, :, :)     ! triplet channel
  
  !single-particle property
  !complex(dp), allocatable :: Gkw(:)           ! Green's function
  complex(dp), allocatable :: Sigma(:)          ! Self-energy
  complex(dp), allocatable :: Sigma_compare(:)     ! Self-energy for comparisons
  
  complex(dp), allocatable :: Sigma_DMFT(:)    ! DMFT Self-energy
  
  complex(dp), allocatable :: chi0_ph(:)       ! bubble diagram in p-h channel
  complex(dp), allocatable :: chi0_pp(:)       ! bubble diagram in p-p channel
  !complex(dp), allocatable :: k1_d(:,:,:)      ! kernel-1 approximation;
  !complex(dp), allocatable :: k1_m(:,:,:)      ! will be used in a future version
  !complex(dp), allocatable :: k1_s(:,:,:)
  !complex(dp), allocatable :: k1_t(:,:,:)
  complex(dp), allocatable :: k2_d1(:,:,:)     ! kernel approximation of Phi_d
  complex(dp), allocatable :: k2_m1(:,:,:)     ! kernel approximation of Phi_m
  complex(dp), allocatable :: k2_s1(:,:,:)     ! kernel approximation of Phi_s
  complex(dp), allocatable :: k2_t1(:,:,:)     ! kernel approximation of Phi_t  
  
  complex(dp), allocatable :: k2_d2(:,:,:)     ! kernel approximation of Phi_d
  complex(dp), allocatable :: k2_m2(:,:,:)     ! kernel approximation of Phi_m
  complex(dp), allocatable :: k2_s2(:,:,:)     ! kernel approximation of Phi_s
  complex(dp), allocatable :: k2_t2(:,:,:)     ! kernel approximation of Phi_t
  
  complex(dp), allocatable :: mat(:,:,:)       ! temperary arry for parquet equation and self-energy calculations
  
  !FFT 
  real(dp), allocatable :: C_wave_x(:)         ! work array for FFT in x
  real(dp), allocatable :: C_wave_y(:)         ! work array for FFT in y
  
  type :: indxmap
    integer :: ix
    integer :: iy
    integer :: iw
  end type indxmap
  
  type(indxmap), allocatable :: Index_fermionic(:), Index_bosonic(:), Index_bosonic_IBZ(:)
  
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
    
    ! ... local vars ...
    integer   :: i, j, k, itemp, idx, ix, ix2, iy, iy2
    real(dp)  :: t1, t2, t3, t4, dkx, dky, kx, ky
      
    real(dp) :: tprime       ! t' in units of t for the t-t' Hubbard model
    
    open(unit=1, file='input.parquet', status='old')
    
    read(1, *) Nx_IBZ, Ny, nOmega, Ngrain, beta, xU, nParticle, tprime
    ! Ny is not needed for the square lattice
    read(1, *) f_damping, f_range, f_sigma
    read(1, *) ite_max  ! maximal number of iterations    
    read(1, *) DGA          ! if the local fully irreducible
                            !vertex from DMFT is used or not
    read(1, *) Sigma_ini    ! if Sigma_ini is provided   
    read(1, *) No_sigma_update  ! For comparisons  
    read(1, *) Eigen    ! If .true. eigenvalues are calculated at every         
                        !iteration (helps to asess stability) 
                    
    close(1)
    !call initialize_lattice(lattice_type) !this will be included later 
    
    Nx = 2*Nx_IBZ-2 ! Lattice dependent
    Ny = 2*Nx_IBZ-2 ! Lattice dependent
    
    Nc = Nx*Ny
    Nf = 2*nOmega
    Nt = Nc*Nf
    Nred = Nx_IBZ*(Nx_IBZ+1)/2*nOmega ! Lattice dependent
    Ns = 8 ! Lattice dependent
    lattice_type = 'Square'
    
    Nb = Nred/ntasks
    
    if (mod(Nred, ntasks) /= 0) then
      call error_msg('The allocated total num. of CPU is not compatible for parallization.')
    else
      if (id == master) then
        write(*, '(a,i4,a)'), 'each core gets ', Nb, ' copies of Nc*Nf data.'
      end if
    end if

    if (Nc == 1) then

    if (.NOT. allocated(Delta)) allocate(Delta(Nf))
    
    open(unit=1, file='Delta.dat', status='unknown')
    do i = Nf/2+1, Nf
      read(1, *) t1, t2, t3
      Delta(i) = dcmplx(t2, t3)
      Delta(Nf-i+1) = conjg(Delta(i))
    end do
    close(1)

    end if
    
    !  --- set up the index map for the convenience of reference ---
    if (.NOT. allocated(Index_fermionic)) allocate(Index_fermionic(Nt))
    if (.NOT. allocated(Index_bosonic))   allocate(Index_bosonic(Nt/2))
    if (.NOT. allocated(Index_bosonic_IBZ))   allocate(Index_bosonic_IBZ(Nred))
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
    
    ! Lattice dependent
    do i = 1, Nx_IBZ
      do j = 1, i
        do k = 1, Nf/2
          idx = (((i-1)*i)/2+j-1)*Nf/2 + k  ! Lattice dependent
          Index_bosonic_IBZ(idx) = indxmap(i, j, k)
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
        if (Ny > 1) Ek(i, j) = Ek(i, j) - Two*cos(ky)*(One+Two*tprime*cos(kx))
      end do
    end do
    
    if (Nc == 1) Ek = Zero
      ! --- set up the tight-binding dispersion for the enlarged case ---
      if (.NOT. allocated(Ek_grain)) allocate(Ek_grain(Nx,Ngrain, Ny,Ngrain))
      dkx = Two*Pi/Nx
      dky = Two*Pi/Ny
      
      Ek_grain = Zero
      do ix = 1, Nx
        do ix2 = 1,Ngrain
          kx = dkx*(ix-1) + (ix2-(Ngrain+1)/2)*dkx/dble(Ngrain)
          do iy = 1, Ny
            do iy2 = 1,Ngrain
              ky = dky*(iy-1)+(iy2-(Ngrain+1)/2)*dky/dble(Ngrain)
              
              !Ek_grain(ix,ix2,iy,iy2) = Ek(ix, iy) 
              if (Nx > 1) Ek_grain(ix,ix2,iy,iy2) = Ek_grain(ix,ix2,iy,iy2)- Two*cos(kx)
              if (Ny > 1) Ek_grain(ix,ix2,iy,iy2) = Ek_grain(ix,ix2,iy,iy2) - Two*cos(ky)*(One+Two*tprime*cos(kx))
              
            end do
          end do            
        end do
      end do  
      
      if (.NOT. allocated(Sigma_H)) allocate(Sigma_H(Nx, Ny))
      Sigma_H = Zero
      
      if (.NOT. allocated(C_wave_x)) allocate(C_wave_x(4*Nx*Ngrain+15))
      if (.NOT. allocated(C_wave_y)) allocate(C_wave_y(4*Ny*Ngrain+15))   
      
      call Zffti(Nx*Ngrain, C_wave_x)
      call Zffti(Ny*Ngrain, C_wave_y)
      
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
      
      if (No_sigma_update) then
        if (.NOT. allocated(Sigma_DMFT))   allocate(Sigma_DMFT(Nf))
        
        open(unit=1, file='Sigma_DMFT.dat', status='old')
        do k = 1, Nf
          read(1, *)  t1, t2, t3
          Sigma_DMFT(k) = dcmplx(t2-Half*xU, t3)        
        end do
        close(1)
      end if        
      
      
   end subroutine readin
  !------------------------------------------------------------------------------
  
    !-------------------------------------------------------
  subroutine pa_Gkw_Chi0_IBZ(ite, Grt)
    !
    ! Purpose
    ! =======
    !   Generate the full Green's function. The bubble diagram of two convoluted Green's function are then 
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
    complex(dp), intent(out) :: Grt(Nx*Ngrain, Ny*Ngrain, Nf)
    
    ! ... local vars ...
    integer     :: i, j, k, idx, iTau, itemp, ix2, iy2, idx1, idx2
    integer     :: i1, j1
    real(dp)    :: t1, t2, t3, mu_UpperBound, mu_LowerBound
    !real(dp)    :: dkx, dky, kx, ky,
    real(dp)    :: w, tau, dummy, FD1, FD2
    complex(dp) :: dummy1D(Nf), coutdata(Nf/2)
    complex(dp) :: Gkw(Nx*Ngrain, Ny*Ngrain, Nf)           ! Green's function
    complex(dp) :: Chi0rt_ph(Nx*Ngrain, Ny*Ngrain, Nf), Chi0rt_pp(Nx*Ngrain, Ny*Ngrain, Nf)
    character(len=30) :: FLE, str1
    character(len=10) :: Mtype
        
    if (.NOT. allocated(Sigma))   allocate(Sigma(Nt))
    if (ite == 1) then 
      
      if (Sigma_ini) then
        open(unit=1, file='Sigma_ini.dat', status='old')
        do i = 1, Nx
          do j = 1, Ny
            do k = 1, Nf
              idx = (Ny*(i-1)+j-1)*Nf + k
              read(1, *) itemp, itemp, t1, t2, t3
              Sigma(idx) = dcmplx(t2, t3)    
            end do
          end do
        end do
        close(1)
      else 
        Sigma = Zero
      end if
      
    end if
    
    if (No_sigma_update) then
      if (.NOT. allocated(Sigma_compare))   allocate(Sigma_compare(Nt))   
      
      Sigma_compare = Sigma
      
      do i = 1, Nx
        do j = 1, Ny
          do k = 1, Nf
            idx = (Ny*(i-1)+j-1)*Nf + k
            Sigma(idx) = Sigma_DMFT(k)
          end do
        end do
      end do    
      
    end if
    
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
    !dkx = Two*Pi/Nx
    !dky = Two*Pi/Ny
    do i = 1, Nx
      do j = 1, Ny
        do ix2 = 1, Ngrain
          do iy2 = 1, Ngrain
            do k = 1, Nf
             w   = Pi/beta*(Two*(k-Nf/2-1)+One)
             idx = (Ny*(i-1)+j-1)*Nf + k
             
             idx1= (i-1)*Ngrain + (ix2-(Ngrain+1)/2) + 1
             if (idx1 < 1) idx1 = idx1 + Ngrain*Nx
             
             idx2= (j-1)*Ngrain + (iy2-(Ngrain+1)/2) + 1
             if (idx2 < 1) idx2 = idx2 + Ngrain*Ny
             
             !if (Nc == 1) then
             ! we switch off single-site functionality
                !Gkw(1,1,k) = One/(xi*w + mu - Delta(k) - Sigma(k))
             !else                
             Gkw(idx1,idx2,k) = One/(xi*w + mu - Ek_grain(i,ix2, j,iy2) - Sigma(idx))
             !end if
             dummy1D(k) = (Gkw(idx1,idx2,k) - One/(xi*w + mu - Ek_grain(i,ix2, j,iy2))) 
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
             Grt(idx1, idx2, iTau) = dummy - &
             exp((beta-tau)*(Ek_grain(i,ix2,j,iy2)-mu))/(One + exp(beta*(Ek_grain(i,ix2, j,iy2)-mu)))
            end do !iTau   
          end do !iy2
        end do !ix2
      end do !j
    end do !i

    ! FFT from k to r space.
    do k = 1, Nf
       call fftf2d(Nx*Ngrain, Ny*Ngrain, Grt(1:Nx*Ngrain, 1:Ny*Ngrain, k), C_wave_x, C_wave_y)
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
          !kx = dkx*(i-1)
          idx1= (i-1)*Ngrain +  1
          do j = 1, Ny
             !  ky = dky*(j-1)
             idx2= (j-1)*Ngrain +  1
             do k = 1, Nf
                write(1, '(2i4, 3f20.12)') i, j, Pi/beta*(Two*(k-Nf/2)-One), Gkw(idx1,idx2,k)
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
    if (.NOT. allocated(Chi0_ph)) allocate(Chi0_ph(Nred))
    if (.NOT. allocated(Chi0_pp)) allocate(Chi0_pp(Nred))
    
    ! --- particle-hole bubble ---
    do k = 1, Nf
      do idx1 = 1, Nx*Ngrain
        if (idx1 == 1) then
          i1 = idx1
        else
          i1 = Nx*Ngrain- idx1 +2
        end if
        
        do idx2 = 1, Ny*Ngrain
          if (idx2 == 1) then
            j1 = 1
          else
            j1 = Ny*Ngrain- idx2 +2
          end if
          Chi0rt_ph(idx1, idx2, k) = -Grt(idx1, idx2, k)*Grt(i1, j1, Nf-k+1)     
        end do
      end do
      call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_ph(1:Nx*Ngrain, 1:Ny*Ngrain, k), C_wave_x, C_wave_y)
    end do
    
    Mtype = 'Bosonic'
    do i = 1, Nx_IBZ
      idx1 = (i-1)*Ngrain +  1
      do j = 1, i
        idx2 = (j-1)*Ngrain +  1
        do k = 1, Nf
          dummy1D(k) = Chi0rt_ph(idx1, idx2, k) 
        end do
        call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
        call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
        do k = 1, Nf/2
          idx = ((i*(i-1))/2+j-1)*Nf/2 + k
          Chi0_ph(idx) = coutdata(k)
        end do
      end do
    end do
    
    ! --- particle-particle bubble ---
    do k = 1, Nf
      do idx1 = 1, Nx*Ngrain
        do idx2 = 1, Ny*Ngrain
          Chi0rt_pp(idx1, idx2, k) = -Half*Grt(idx1, idx2, k)**2
        end do
      end do
      call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_pp(1:Nx*Ngrain, 1:Ny*Ngrain, k), C_wave_x, C_wave_y)
    end do
    
    do i = 1, Nx_IBZ
      idx1 = (i-1)*Ngrain + 1
      do j = 1, i
        idx2 = (j-1)*Ngrain +  1
        do k = 1, Nf
          dummy1D(k) = Chi0rt_pp(idx1, idx2, k)    
        end do
        
        call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
        call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
        
        do k = 1, Nf/2
          idx = ((i*(i-1))/2+j-1)*Nf/2 + k
          Chi0_pp(idx) = coutdata(k)
        end do          
      end do
    end do
    
    !--- monitoring the bubble results ---
    if (id == master) then
      FLE = 'Chi0-'//trim(str1)//'.dat'
      open(unit=1, file=FLE, status='unknown')
      do i = 1, Nx_IBZ
        !kx = dkx*(i-1)
        do j = 1, i
          !ky = dky*(j-1)
          do k = 1, Nf/2
            idx = ((i*(i-1))/2+j-1)*Nf/2 + k
            write(1, '(2i4, 5f12.6)') i, j, Pi/beta*Two*(k-1), &
                  Chi0_ph(idx), Chi0_pp(idx) 
          end do
          write(1, *)
        end do
      end do
      close(1)
    end if
    
  end subroutine pa_Gkw_Chi0_IBZ
  !------------------------------------------------------------------------------  
  
  !------------------------------------------------------------------------------
  subroutine list_symmetries(typ, map, symm_list)
    
    !---------Purpose-------------
    !  This subroutine lists all operations that need to be performed to create
    !  star k for a given k. We need it to restore full BZ from IBZ
    !  The operations are numbered as follows (now square lattice only):
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx)
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx)
    
    character(len=30), intent(in) :: typ
    type(indxmap), intent(in) :: map
    logical, intent(out) :: symm_list(Ns) 
    
    integer ::  kx, ky 
    
    kx=map%ix
    ky=map%iy
    
    symm_list=.False.
    
    if (((kx==1).and.(ky==1)).or.((kx==Nx_IBZ).and.(ky==Nx_IBZ))) then ! (0,0),(Pi,Pi)
      
      symm_list(1)=.True. ! identity
      
    else if  ((kx==Nx_IBZ).and.(ky==1)) then  !(Pi,0) (2x)
      
      symm_list(1)=.True. ! identity
      symm_list(5)=.True. ! 5  (kx,ky) -> (ky,kx)
      
    else if  ((kx==Nx_IBZ).and.((ky>1).and.(ky<Nx_IBZ))) then !(Pi,j<Pi) (4x)
      
      symm_list(1)=.True. ! identity
      symm_list(2)=.True.  !  2  (kx,ky) -> (kx,-ky)
      symm_list(5)=.True. ! 5  (kx,ky) -> (ky,kx)
      symm_list(8)=.True.  !  8  (kx,ky) -> (-ky, kx)
      
    else if ((kx>1).and.(ky<Nx_IBZ)) then
      if (ky==1) then !(0<i<Pi,0) (4x)
        
        symm_list(1)=.True. ! identity
        symm_list(3)=.True.  !  3  (kx,ky) -> (-kx,ky)
        symm_list(5)=.True.  !  5  (kx,ky) -> (ky,kx)
        symm_list(6)=.True.  !  6  (kx,ky) -> (ky, -kx)
        
      else if (ky==kx) then ! (diagonal) (4x)
        
        symm_list(1)=.True. ! identity
        symm_list(2)=.True.  !  2  (kx,ky) -> (kx,-ky) 
        symm_list(3)=.True.  !  3  (kx,ky) -> (-kx,ky)   
        symm_list(4)=.True.  !  4  (kx,ky) -> (-kx,-ky)
        
      else ! all other points (8x)
        
        symm_list=.True.
        
      end if
    end if
    
  end subroutine list_symmetries
  !------------------------------------------------------------------------------
  
  subroutine check_sym(i_in, i_sym, i_out)

    integer, intent(out) :: i_sym
    type(indxmap), intent(in) :: i_in
    type(indxmap), intent(out) :: i_out
    
    integer kx, ky, kw, jx, jy
    
    !---------Purpose-------------
    !  This subroutine checks inverse of which symmetry operation needs to be performed
    !  to bring the momentum from full BZ into IBZ and performs this inverse operation on i_in
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx) inv of 8
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx) inv of 6

    kx = i_in%ix
    ky = i_in%iy
    kw = i_in%iw
    
    if (kx >= ky) then
        if (kx <= Nx_IBZ) then 
            i_sym = 1 
            i_out = i_in            
        else if (ky >= Nx_IBZ) then
            i_sym = 7
            jx= -ky + Ny + 2
            if (jx > Nx) jx = jx - Nx
            jy= -kx + Nx + 2
            if (jy > Ny) jy = jy - Ny
            i_out=indxmap(jx,jy,kw)    
        else if (ky >= Nx + 2 - kx) then
            i_sym = 8
            jx = ky
            jy= -kx + Nx + 2
            if (jy > Ny) jy = jy - Ny
            i_out=indxmap(jx,jy,kw)        
        else 
            i_sym = 3
            jx = -kx + Nx + 2
            if (jx > Nx) jx = jx - Nx
            jy = ky
            i_out=indxmap(jx,jy,kw)
        end if                
    else if (kx < ky) then                        
        if (kx >= Nx_IBZ) then 
            i_sym = 4
            jx = -kx + Nx + 2
            if (jx > Nx) jx = jx - Nx
            jy = -ky + Ny + 2
            if (jy > Ny) jy = jy - Ny
            i_out=indxmap(jx,jy,kw)        
        else if (ky <= Nx_IBZ) then 
            i_sym = 5
            jx = ky
            jy = kx
            i_out=indxmap(jx,jy,kw)
        else if (kx >= Ny + 2 - ky) then
            i_sym = 2        
            jx = kx
            jy = -ky + Ny + 2
            if (jy > Ny) jy = jy - Ny
            i_out=indxmap(jx,jy,kw)            
        else 
            i_sym = 6     
            jx= -ky + Ny + 2
            if (jx > Nx) jx = jx - Nx
            jy = kx
            i_out=indxmap(jx,jy,kw)
        end if                       
    end if      
    
end subroutine check_sym
  
  
  !------------------------------------------------------------------------------
  subroutine symmetry_operation(typ, i_sym, i_in, i_out)
    !
    !  Symmetry operation on one index map.
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx) inv of 8
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx) inv of 6
    
    character(len=30), intent(in) :: typ
    integer, intent(in) :: i_sym
    type(indxmap), intent(in) :: i_in
    type(indxmap), intent(out) :: i_out
    
    integer :: ix, iy, wi, kx, ky
    
    ix=i_in%ix
    iy=i_in%iy
    wi=i_in%iw 
    
    select case (i_sym)
      case (1)
        i_out=i_in
      !if (typ == 'Square') then
      case (2)  ! (kx,ky) -> (kx,-ky)
        kx = ix
        ky = -iy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
      case (3) ! (kx,ky) -> (-kx,ky)
        kx = -ix + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = iy
        i_out=indxmap(kx,ky,wi)
      case (4)  ! (kx,ky) -> (-kx,-ky)
        kx = -ix + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = -iy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
      case (5) ! (kx,ky) -> (ky,kx)
        kx=iy
        ky=ix
        i_out=indxmap(kx,ky,wi)
      case (6)  !  (kx,ky) -> (ky,-kx)
        kx=iy
        ky= -ix + Nx + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
      case (7)  ! (kx,ky) -> (-ky,-kx)
        kx= -iy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= -ix + Nx + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
      case (8)  ! (kx,ky) -> (-ky, kx)
        kx= -iy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= ix
        i_out=indxmap(kx,ky,wi)
    end select
    !end if
    
  end subroutine symmetry_operation
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine symmetry_operation_inv(typ, i_sym, i_in, j_in, i_out, j_out)
    !
    !  Inverse symmetry operation on two index maps.
    !
    !  1  identity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx) inv of 8
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx) inv of 6
    
    character(len=30), intent(in) :: typ
    integer, intent(in) :: i_sym
    type(indxmap), intent(in) :: i_in, j_in 
    type(indxmap), intent(out) :: i_out, j_out
    
    integer :: ix, iy, jx, jy, wi, wj, kx, ky
    
    ix=i_in%ix
    iy=i_in%iy
    wi=i_in%iw 
    
    jx=j_in%ix
    jy=j_in%iy
    wj=j_in%iw  
    
    
    
    select case (i_sym)
      case (1)
        i_out=i_in
        j_out=j_in
        ! if (typ == 'Square') then
      case (2)  ! (kx,ky) -> (kx,-ky)
        kx = ix
        ky = -iy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
        kx = jx
        ky = -jy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        j_out=indxmap(kx,ky,wj)
      case (3) ! (kx,ky) -> (-kx,ky)
        kx = -ix + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = iy
        i_out=indxmap(kx,ky,wi)
        kx = -jx + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = jy
        j_out=indxmap(kx,ky,wj)
      case (4)  ! (kx,ky) -> (-kx,-ky)
        kx = -ix + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = -iy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
        kx = -jx + Nx + 2
        if (kx > Nx) kx = kx - Nx
        ky = -jy + Ny + 2
        if (ky > Ny) ky = ky - Ny
        j_out=indxmap(kx,ky,wj)
      case (5) ! (kx,ky) -> (ky,kx)
        kx=iy
        ky=ix
        i_out=indxmap(kx,ky,wi)
        kx=jy
        ky=jx
        j_out=indxmap(kx,ky,wj)
      case (6)  ! inverse of 6 is 8 (kx,ky) -> (-ky, kx)
        kx= -iy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= ix
        i_out=indxmap(kx,ky,wi)
        kx= -jy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= jx
        j_out=indxmap(kx,ky,wj)
      case (7)  ! (kx,ky) -> (-ky,-kx)
        kx= -iy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= -ix + Nx + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
        kx= -jy + Ny + 2
        if (kx > Nx) kx = kx - Nx
        ky= -jx + Nx + 2
        if (ky > Ny) ky = ky - Ny
        j_out=indxmap(kx,ky,wj)
      case (8)  ! inverse of 8 is 6: (kx,ky) -> (ky,-kx)
        kx=iy
        ky= -ix + Nx + 2
        if (ky > Ny) ky = ky - Ny
        i_out=indxmap(kx,ky,wi)
        kx=jy
        ky= -jx + Nx + 2
        if (ky > Ny) ky = ky - Ny
        j_out=indxmap(kx,ky,wj)
    end select
    
    !end if
    
  end subroutine symmetry_operation_inv
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  integer function list_index_F(P)  result(idx)
    type(Indxmap), intent(in) :: P
    idx =((P%ix-1)*Ny+P%iy-1)*Nf+P%iw
  end function list_index_F
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  integer function list_index(P, typ) result(idx)
  ! this function is obsolete in version 1.1
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
  
  !------------------------------------------------------------------------------
  integer function list_index_IBZ(P, typ) result(idx)
    type(Indxmap), intent(in) :: P
    character(len=30), intent(in) :: typ
    
    !(only bosonic is now in the IBZ)
    idx = (((P%ix-1)*P%ix)/2+P%iy-1)*Nf/2+P%iw 
    
  end function list_index_IBZ
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  integer function list_index_B(P) result(idx)
    !now the same as list_index_IBZ, one of them will be made redundant soon
    type(Indxmap), intent(in) :: P
    idx = (((P%ix-1)*P%ix)/2+P%iy-1)*Nf/2+P%iw 
  end function list_index_B
  !------------------------------------------------------------------------------
  
  
  !-------------------------------------------------------
  subroutine index_FaddB(idx1, idx2, final_Indx)
    !
    ! Purpose
    ! =======
    !  For given two index in the complete list, find out the resulting index for
    !  given operation.
    !
    type(Indxmap), intent(in)     :: idx1, idx2
    type(Indxmap), intent(out)    :: final_Indx
    
    ! ... local vars ...
    integer :: i, j, k
    
    i = idx1%ix + idx2%ix - 1
    if (i > Nx) i = i - Nx
    j = idx1%iy + idx2%iy - 1
    if (j > Ny) j = j - Ny
    k = idx1%iw + idx2%iw - 1      
    final_Indx = indxmap(i, j, k)
   
  end subroutine index_FaddB
  !-------------------------------------------------------
  
  !-------------------------------------------------------
  subroutine index_FaddF(idx1, idx2, final_Indx)
    !
    ! Purpose
    ! =======
    !  For given two index in the complete list, find out the resulting index for
    !  given operation.
    !
    type(Indxmap), intent(in)     :: idx1, idx2
    type(Indxmap), intent(out)    :: final_Indx
    
    ! ... local vars ...
    integer :: i, j, k
    
  
    i = idx1%ix + idx2%ix - 1
    if (i > Nx) i = i - Nx
    j = idx1%iy + idx2%iy - 1
    if (j > Ny) j = j - Ny
    k = idx1%iw + idx2%iw - Nf 
    final_Indx = indxmap(i, j, k)    

   
  end subroutine index_FaddF
  !------------------------------------------------------- 

  !-------------------------------------------------------
  subroutine index_minusF(idx1, final_Indx)
    !
    ! Purpose
    ! =======
    !  For given two index in the complete list, find out the resulting index for
    !  given operation.
    !
    type(Indxmap), intent(in)     :: idx1
    type(Indxmap), intent(out)    :: final_Indx
    
    ! ... local vars ...
    integer :: i, j, k
    
    i = -idx1%ix + Nx + 2
    if (i > Nx) i = i - Nx
    j = -idx1%iy + Ny + 2
    if (j > Ny) j = j - Ny
    k = -idx1%iw + Nf + 1
    final_Indx = indxmap(i, j, k)
    
  end subroutine index_minusF
  !-------------------------------------------------------
    !-------------------------------------------------------
  subroutine index_minusB(idx1, final_Indx)
    !
    ! Purpose
    ! =======
    !  For given two index in the complete list, find out the resulting index for
    !  given operation.
    !
    type(Indxmap), intent(in)     :: idx1
    type(Indxmap), intent(out)    :: final_Indx
    
    ! ... local vars ...
    integer :: i, j, k
    
    i = -idx1%ix + Nx + 2
    if (i > Nx) i = i - Nx
    j = -idx1%iy + Ny + 2
    if (j > Ny) j = j - Ny
    k = -idx1%iw + 2
    final_indx = indxmap(i, j, k)
    
  end subroutine index_minusB
  !-------------------------------------------------------

  !------------------------------------------------------------------------------
  subroutine output_2k(data1, data2, k, name, unt)
    !
    ! Purpose
    ! =======
    !   Write data with two four-k dependencies to text file 
    !   specified by "name" with unit number "unit"
    !
    complex(dp), allocatable , intent(in) :: data1(:, :, :),data2(:, :, :)
    character(len = *), intent(in)        :: name
    integer,     intent(in)               :: k, unt          
    integer :: ix2, iy2, iw2, ix, iy, iw, idx2, idx
    
    open(unit=unt, file=trim(name), status='unknown')
    
    do ix2 = 1,Nx
      do iy2= 1,Ny
        do iw2 = 1,Nf
          idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
          do ix = 1,Nx
            do iy= 1,Ny
              do iw = 1,Nf
                idx = ((ix-1)*Ny+iy-1)*Nf + iw       
                
                write(unt, '( 6i5, 4f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(data1(idx2,idx,k)), aimag(data1(idx2,idx,k)), dble(data2(idx2,idx,k)), aimag(data2(idx2,idx,k))
              end do
            end do
          end do
        end do
      end do
    end do             
    
    close(unt) 
    
  end subroutine output_2k
  
    !------------------------------------------------------------------------------
  subroutine output_2k_1(data1, k, name, unt)
    !
    ! Purpose
    ! =======
    !   Write data with two four-k dependencies to text file 
    !   specified by "name" with unit number "unit"
    !
    complex(dp), allocatable , intent(in) :: data1(:, :, :)
    character(len = *), intent(in)        :: name
    integer,     intent(in)               :: k, unt          
    integer :: ix2, iy2, iw2, ix, iy, iw, idx2, idx
    
    open(unit=unt, file=trim(name), status='unknown')
    
    do ix2 = 1,Nx
      do iy2= 1,Ny
        do iw2 = 1,Nf
          idx2 = ((ix2-1)*Ny+iy2-1)*Nf + iw2
          do ix = 1,Nx
            do iy= 1,Ny
              do iw = 1,Nf
                idx = ((ix-1)*Ny+iy-1)*Nf + iw       
                
                write(unt, '( 6i5, 2f15.8)') ix2, iy2, iw2-Nf/2-1,ix, iy, iw-Nf/2-1, dble(data1(idx2,idx,k)), aimag(data1(idx2,idx,k))
              end do
            end do
          end do
        end do
      end do
    end do             
    
    close(unt) 
    
  end subroutine output_2k_1

  !------------------------------------------------------------------------------
  subroutine Memory_Release
    
    if (allocated(Delta)) Deallocate(Delta)
    
    ! clear the dispersion 
    if (allocated(Ek)) Deallocate(Ek)
    if (allocated(Ek_grain)) Deallocate(Ek_grain)
    
    ! deallocate the Hatree energy
    if (allocated(Sigma_H)) Deallocate(Sigma_H)
    
    ! deallocate the index array
    if (allocated(Index_bosonic)) Deallocate(Index_bosonic)
    if (allocated(Index_bosonic_IBZ)) Deallocate(Index_bosonic_IBZ)
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
    ! if (allocated(Gkw))     Deallocate(Gkw)
    if (allocated(Sigma))   Deallocate(Sigma)
    !if (allocated(Gloc))   Deallocate(Gloc)
    if (allocated(Sigma_DMFT))   Deallocate(Sigma_DMFT)
    
    if (allocated(Sigma_compare))   Deallocate(Sigma_compare)   
    if (allocated(Chi0_ph)) Deallocate(Chi0_ph)
    if (allocated(Chi0_pp)) Deallocate(Chi0_pp)
    
    !if (allocated(K1_d))   Deallocate(K1_d)
    !if (allocated(K1_m))   Deallocate(K1_m)
    !if (allocated(K1_s))   Deallocate(K1_s)
    !if (allocated(K1_t))   Deallocate(K1_t)
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
  !------------------------------------------------------------------------------  
  
end module parquet_util
