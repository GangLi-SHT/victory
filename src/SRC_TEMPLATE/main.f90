program main
  
  use mpi
  use global_parameter

  ! ... local vars ...
  integer :: ite
 
  ! ------------------- initialize the mpi enviroment ----------------------
  call parallel_start

  ! ------------------- print the license message --------------------------
  master = 0
  if (id ==master) call license_message(ntasks)

  ! ------------------- start the self-consisent loop ----------------------
  do ite = 1, 10
     if (id == master) call loop_message(ite)

     
  end do

  ! ------------------- finalize the mpi enviroment ------------------------
  call parallel_end

end program main
