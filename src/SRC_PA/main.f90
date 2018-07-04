program main
  
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util
  use parquet_kernel
  use parquet_equation
  use parquet_selfenergy
  use parquet_phys
  
  implicit none
  
  ! ... local vars ...
  integer       :: ite
  real(dp)      :: t1, t2, t3, t4, t5, t6
  logical       :: Converged
  
  complex(dp), allocatable :: Grt(:, :, :)
  
  ! ------------------- initialize the mpi enviroment ----------------------
  call parallel_start
  
  ! ------------------- print the license message --------------------------
  master = 0
  if (id == master) call license_message(ntasks)
  
  call readin
  
  allocate(Grt(Nx, Ny, Nf))
  ! ------------------- start the self-consisent loop ----------------------
  Converged = .False.
  ite = 1
  do while (.NOT. Converged .and. ite < ite_max)
    
    if (id == master) call loop_message(ite)
    call cpu_time(t1)
    
    ! calculate single-particle green's function
    call pa_Gkw_Chi0_IBZ(ite, Grt)
    call cpu_time(t2)
    
    if (id == master) write(*, "(a, f12.6)") '  time spent on pa_Gkw_Chi0 is:', t2-t1 
    
    ! determine reducible vertex function and its kernel approximation
    !if (edge_scan) then
      call reducible_vertex(ite)
      call get_kernel_function(ite)
    !else
      !call get_kernel_function(ite)
      !call reducible_vertex(ite)
    !end if 
    call cpu_time(t3)
    if (id == master) write(*, "(a, f12.6)") '  time spent on kernel calculation is:', t3-t2
    
    ! solve the parquet equation
    call solve_parquet_equation(ite)
    call cpu_time(t4)
    if (id == master) write(*, "(a, f12.6)") '  time spent on solve_parquet_equation is:', t4-t3
    
    ! calculate the self-energy 
    call self_energy(ite, Grt, converged)
    call cpu_time(t5)
    if (id == master) write(*, "(a, f12.6)") '  time spent on self_energy is:', t5-t4
    
    ! --- update irreducible vertex in each channel ---
    G_d = F_d - G_d   
    G_m = F_m - G_m
    G_s = F_s - G_s
    G_t = F_t - G_t
    
    ! ------------------- calculate eigen values in each channel -------------------------------
        
    call solve_eigen_equation(ite,converged)
    
    call cpu_time(t6) 
    if (id == master) write(*, "(a, f12.6)") '  time spent on eigen_equation is:', t6-t5
    
    if (id == master) write(*, "(a, f12.6)") 'total time cost:', t6-t1
    ite = ite + 1
  end do
  
  
  ! ------------------- clean the memory -----------------------------------
  !call MPI_barrier(MPI_COMM_WORLD, rc) !not needed
  call Memory_Release
  
  ! ------------------- finalize the mpi enviroment ------------------------
  call parallel_end
  
  deallocate(Grt)
end program main
