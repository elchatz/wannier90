module aux
  integer, parameter :: dp = kind(1.0d0)
  private dp
contains
  subroutine read_eval(nb, nk, eval)
    integer, intent(in) :: nb, nk
    real(kind=dp), allocatable, intent(inout) :: eval(:, :)
    allocate (eval(nb, nk))
    open (newunit=iu, file='gaas.eig')
    do ik = 1, nk
    do ib = 1, nb
      read (iu, *) ijunk, ijunk, eval(ib, ik)
    enddo
    enddo
  end subroutine read_eval

  subroutine read_kp(nk, kpcart)
    integer, intent(in) :: nk
    real(kind=dp), allocatable, intent(out) :: kpcart(:, :)
    allocate (kpcart(3, nk))
    open (newunit=iu, file='KPOINTS')
    do ik = 1, nk
      read (iu, *) (kpcart(ic, ik), ic=1, 3)
    enddo
  end subroutine read_kp

  subroutine read_cell(nbas, uccart, xcart)
    character(len=100) :: junk
    integer, intent(inout) :: nbas
    real(kind=dp), allocatable, intent(inout) :: xcart(:, :)
    real(kind=dp), intent(inout) :: uccart(:, :)
    open (newunit=iu, file='CELL')
    do ib = 1, 3
      read (iu, *) (uccart(ic, ib), ic=1, 3)
    enddo
    open (newunit=iu, file='POSITIONS')
    read (iu, *) nbas
    allocate (xcart(3, nbas))
    do ib = 1, nbas
      read (iu, '(a)') junk
    enddo
    do ib = 1, nbas
      read (iu, *) (xcart(ic, ib), ic=1, 3)
    enddo
  end subroutine read_cell
end module

program ok
  use aux
  use mpi_f08
  use w90_helper_types
  use w90_comms, only: w90_comm_type

  implicit none

  complex(kind=dp), allocatable :: a(:, :, :)
  complex(kind=dp), allocatable :: m(:, :, :, :)
  complex(kind=dp), allocatable :: morig(:, :, :, :)
  complex(kind=dp), allocatable :: u(:, :, :)
  complex(kind=dp), allocatable :: uopt(:, :, :)
  integer :: nbas, ierr, stdout, stderr
  integer :: nb, nk, nn, nw, nkabc(3), exclude(5)
  integer, allocatable :: distk(:)
  real(kind=dp), allocatable :: eval(:, :)
  real(kind=dp), allocatable :: kpcart(:, :) ! cartesian k list
  real(kind=dp), allocatable :: xcart(:, :) ! cartesian atom positions
  real(kind=dp) :: uccart(3, 3) ! cartesian unit cell

  type(lib_global_type), target :: w90main
  type(lib_w90_type), target :: w90dat
  type(w90_comm_type) :: comm

  nb = 12
  nw = 8
  nkabc = (/4, 4, 4/)
  exclude = (/1, 2, 3, 4, 5/)
  nk = nkabc(1)*nkabc(2)*nkabc(3)

  ! MPI
  comm%comm = mpi_comm_world
  call mpi_init(ierr)

  ! io and k distribution
  call get_fortran_stdout(stdout)
  call get_fortran_stderr(stderr)
  allocate (distk(nk))
  distk = 0
  call set_kpoint_distribution(w90main, distk, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  ! required settings
  call set_option(w90main, 'num_wann', nw)
  call set_option(w90main, 'num_bands', nb)
  call set_option(w90main, 'num_kpts', nk)
  call read_cell(nbas, uccart, xcart)
  call set_option(w90main, 'unit_cell_cart', uccart)
  call set_option(w90main, 'mp_grid', nkabc)
  call read_kp(nk, kpcart)
  call set_option(w90main, 'kpoints', kpcart)
  !call set_option(w90main, 'dis_froz_max', 14.0d0) ! not optional for disentanglement? fixme

  call set_io_seedname(w90main, 'gaas')

  ! apply settings (and discard settings store)
  call input_setopt(w90main, w90dat, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  ! optional settings
  ! collect other settings from hacked win file
  call set_io_seedname(w90main, 'gaas-edit')
  call input_reader(w90main, w90dat, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  ! go back to original extension
  call set_io_seedname(w90main, 'gaas')

  ! k setup needs attention
  ! must be done before reading overlaps
  call create_kmesh(w90main, stdout, stderr, ierr, comm)
  nn = w90main%kmesh_info%nntot
  write (*, *) "nb, nw, nk, nn, np: ", nb, nw, nk, nn, w90dat%num_proj

  ! setup all matrices
  allocate (a(nb, nw, nk))
  call set_a_matrix(w90dat, a)
  allocate (morig(nb, nb, nn, nk))
  call set_m_orig(w90dat, morig) ! m_matrix_local_orig
  allocate (m(nw, nw, nn, nk)) ! m also known as m_matrix_local
  call set_m_matrix_local(w90dat, m)
  allocate (u(nw, nw, nk))
  call set_u_matrix(w90main, u)
  allocate (uopt(nb, nw, nk))
  call set_u_opt(w90main, uopt)

  ! read from ".mmn" and ".amn"
  ! and assign to m_orig and a
  ! (or m and u if not disentangling)
  call overlaps(w90main, w90dat, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  call read_eval(nb, nk, eval)
  call set_eigval(w90main, eval)

  call disentangle(w90main, w90dat, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  ! change a setting
  call set_option(w90main, 'iprint', 2)
  call input_setopt(w90main, w90dat, stdout, stderr, ierr, comm)

  call wannierise(w90main, w90dat, stdout, stderr, ierr, comm)
  if (ierr /= 0) stop

  call mpi_finalize()
end program
