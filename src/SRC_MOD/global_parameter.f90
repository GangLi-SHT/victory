module global_parameter
  
  implicit none

  integer,  parameter :: dp = selected_real_kind(8)
  real(dp), parameter :: pi = acos(-1.d0)
  real(dp), parameter :: Zero = 0.d0, Half = 0.5d0, One = 1.d0, Two = 2.d0  
  complex(dp), parameter :: xi = dcmplx(Zero, One)
  complex(dp), parameter :: Zero_c = dcmplx(Zero, Zero), One_c = dcmplx(One, Zero)
  complex(dp), parameter :: Half_c = dcmplx(Half, Zero), Two_c  = dcmplx(Two, Zero)

  !--- mpi ---
  integer  :: ntasks, master, id, rc
  integer, allocatable :: status(:), send_request(:), recv_request(:)

  integer, parameter :: nSpin   = 2

  integer :: nOmega  = 64 
  integer :: nTau    = 200

  real(dp) :: beta = 2.d0
  real(dp) :: xU   = 4.d0
  real(dp) :: xJ   = 1.d0
  real(dp) :: mu   = 0.d0
  real(dp) :: nParticle = 1.d0
  
  complex(dp), allocatable :: Omega(:)

  logical :: Debug_Info = .False.

     ! compute an LU factorization of a general M-by-N matrix A using
     ! partial pivoting with row interchanges
  interface
     subroutine ZGETRF( M, N, A, LDA, IPIV, INFO )
       integer, intent(in)    :: LDA, M, N
       integer, intent(out)   :: IPIV(*), INFO
       complex(Selected_Real_Kind(8)), intent(inout) :: A(LDA, N)
     end subroutine ZGETRF

     ! compute  the  inverse of a matrix using the LU factorization
     ! computed by ZGETRF

     subroutine ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       integer, intent(in)    :: N, LDA, LWORK
       integer, intent(in)    :: IPIV(*)
       integer, intent(out)   :: INFO
       complex(Selected_Real_Kind(8)), intent(inout)  :: A(LDA, N)
       complex(Selected_Real_Kind(8)), intent(out)    :: WORK(LWORK)
     end subroutine ZGETRI

     subroutine DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
       character, intent(in) :: JOBZ
       integer, intent(in)   :: LDZ, N
       integer, intent(out)  :: INFO
       real(selected_real_kind(8)), intent(inout) :: D( * ), E( * )
       real(selected_real_kind(8)), intent(out)   :: WORK( * ), Z( LDZ, * )
     end subroutine DSTEV

     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

       INTEGER, intent(in)   ::    INFO, LDA, LDB, N, NRHS
       INTEGER, intent(out)  ::    IPIV( * )
       real(selected_real_kind(8)), intent(in)  :: A( LDA, * )
       real(selected_real_kind(8)), intent(out) :: B( LDB, * )
     END SUBROUTINE DGESV
    end interface

    interface inverse
       module procedure inverse_cmplx, inverse_dble
    end interface inverse

    interface get_userARG
       module procedure get_userARG_int, get_userARG_dble
    end interface get_userARG

contains
  !----------------------------------------------------------------------------------------------
  subroutine license_message(num_node)
    implicit none
    integer, intent(in) :: num_node

    write(*, '(a)') ''//achar(27)//'[33m +--------------------------------------------------+'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m         _                                      '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m        (_)      _                              '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m   _   _ _  ____| |_  ___   ____ _   _          '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m  | | | | |/ ___)  _)/ _ \ / ___) | | |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m   \ V /| ( (___| |_| |_| | |   | |_| |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m    \_/ |_|\____)\___)___/|_|    \__  |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m                                (____/          '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m    Vienna Computational Tool Depository        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[34m    ==     =             ==           ==        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  If you find victory useful in your research,  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   you should                                   '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   (i)  cite our papers and the appropriate     '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        references therein AND                  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   (ii) state in your manuscript/paper that you '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        have used the victory code (or a        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        modified version of it). An appropriate '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        way of acknowledging the use of victory '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        in your publications would be, e.g.     '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        adding a sentence like                  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m "The calculation has been performed using the  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  victory code", followed by the citation to our'//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  paper: Phys. Rev. B 93, 165103 (2016).        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m             Gang Li, Karsten Held @ copyright  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                    2015.06.01  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m +--------------------------------------------------+'//achar(27)//'[0m'
                       
    write(*, *)        
    write(*, '(2x, a, i5, a)') '... MPI starts with', num_node,'  nodes. ...'

  end subroutine license_message
  
  !----------------------------------------------------------------------------------------------  
  subroutine loop_message(ite)
    implicit none
    integer, intent(in) :: ite
    
    ! ---- self-consistent loop monitoring ---
    write(*,*)
    write(*,"(a36)") '!----------------------------------!'
    write(*,"(a24,i5,a7)") '!   Self-Consistent loop', ite, '      !'
    write(*, "(a36)") '!----------------------------------!'
    
  end subroutine loop_message

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine error_msg(msg)

    character(len=*), intent(in) :: msg

    write(*,*) trim(msg)
    stop
  end subroutine error_msg

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine get_userARG_int(arg_message, arg_value)
    !
    ! get user defined values interactively given from, typically, command line. This is used extensively in every algorithm contained
    ! in this pacakge to specify interesting parameters for different problems. If the given parameter is reasonable and accepted, the 
    ! default value will be replaced. Otherwise, it keeps the default one.

    implicit none
    character(len=*), intent(in) :: arg_message
    integer, intent(out) :: arg_value

    character(len=80) :: buffer

    arg_value = 0
    write(*, *)
102 write(*, 101) trim(arg_message)

    read(*, *) arg_value
 
    write(*, '(a15, i8, a35)') 'your inputs are', arg_value, ', do you want to continue? (yes/no)'
    read(*, *) buffer
    if (trim(buffer)=='n' .or. trim(buffer)=='no' .or. trim(buffer)=='N' .or. trim(buffer)=='NO') goto 102 
   
101 format (a<len(arg_message)>)

  end subroutine  get_userARG_int

  subroutine get_userARG_dble(arg_message, arg_value)

    implicit none
    character(len=*), intent(in) :: arg_message
    real(dp), intent(out) :: arg_value

    character(len=80) :: buffer

    arg_value = Zero
    write(*, *)
102 write(*, 101) trim(arg_message)

    read(*, *) arg_value

    write(*, '(a23, f10.6, a35)') 'the input values are', arg_value, ', do you want to continue? (yes/no)'
    read(*, *) buffer
    if (trim(buffer)=='n' .or. trim(buffer)=='no' .or. trim(buffer)=='N' .or. trim(buffer)=='NO') goto 102

101 format (a<len(arg_message)>)

  end subroutine  get_userARG_dble
  
  !------------------------------------------------------------------------------------------------------ 
  subroutine system_mem_usage
    USE IFPORT
    implicit none
    real(kind=8) :: valueRSS

    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer :: pid, itemp, i
    logical :: ifxst
    
    valueRSS=-1    ! return negative number if not found
    
    !--- get process ID
    
    pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'
    
    !--- read system file
    
    inquire (file=filename,exist=ifxst)
    if (.not.ifxst) then
       write (*,*) 'system file does not exist'
       return
    endif
    
    open(unit=100, file=filename, action='read')
    do
       read (100,'(a)',end=120) line
       if (line(1:6).eq.'VmRSS:') then
          read (line(7:),*) itemp
          exit
       endif
    enddo
120 continue

    valueRSS = dble(itemp)
    i = 0
    do while (Int(itemp/1024**i) > 100)
       i = i + 1
    end do
    valueRSS = valueRSS/dble(1024**i)    

    select case (i)
    case (0)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'KB'
    case (1)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'MB'
    case (2)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'GB'
    case (3)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'TB'
    case default
      write(*, '(a)') 'Memory used more than 1024 TB'
    end select
    close(100)

  end subroutine system_mem_usage

  !---------------------------------------------------------------------------------------------------------------------------------
  real(dp) function ranw()
    !
    ! Purpose
    ! =======
    !   random number generator, use the f90 intrinsic routin, random numbers
    !   distuributes in (0, 1)

    call random_number(ranw)

  end function ranw

  !-------------------------------------------------------------------------------------------------------------------
  subroutine inverse_Cmplx(matrix_a, matrix_b, n)
    !
    ! Purpose
    ! =======
    !   complex matrix inversion operation
    !
    ! Arguments
    ! =========
    !
    implicit none
    integer, intent(in) :: n
    complex(dp), intent(in)  :: matrix_a(n, n)
    complex(dp), intent(out) :: matrix_b(n, n)

    integer :: ipiv(n), info
    complex(dp) :: y(n, n), work(n)

    matrix_b = Zero
    y = matrix_a
    call zgetrf(n,n,y,n,ipiv,info)
    call zgetri(n,y,n,ipiv,work,n,info)
    matrix_b = y

  end subroutine inverse_Cmplx

  !--------------------------------------------------------------------------------------------------------------------
  subroutine inverse_dble(matrix_a, matrix_b, n)
    !
    ! Purpose
    ! =======
    !   complex matrix inversion operation
    !
    ! Arguments
    ! =========
    !
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in)  :: matrix_a(n, n)
    real(dp), intent(out) :: matrix_b(n, n)

    integer :: ipiv(n), info
    real(dp) :: y(n, n), work(n)

    matrix_b = Zero
    y = matrix_a
    call dgetrf(n,n,y,n,ipiv,info)
    call dgetri(n,y,n,ipiv,work,n,info)
    matrix_b = y
  end subroutine inverse_dble

end module global_parameter
