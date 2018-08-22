!######################################################################
! FGH1D: A simple program to solve the 1D time-independent Schrodinger
!        equation using the Fourier Grid Hamiltonian method
!######################################################################

module global

  implicit none

  integer, parameter :: dp=selected_real_kind(8)

  ! Number of grid points
  integer :: Npts

  ! Grid bounds
  real(dp), dimension(2) :: bounds

  ! Grid spacing
  real(dp) :: dx
  
  ! Name of the model Hamiltonian
  character(len=20) :: aham

  ! Model Hamiltonian number
  integer :: iham

  ! Hamiltonian parameters
  integer               :: nvpar
  real(dp), allocatable :: vpar(:)
  real(dp)              :: mass
  
  ! Hamiltonian matrix
  real(dp), allocatable :: hmat(:,:)

  ! Eigenvectors and eigenvalues
  real(dp), allocatable :: eigvec(:,:)
  real(dp), allocatable :: eigval(:)
  
end module global

!######################################################################

program fgh1d

  use global
  
  implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  call initialise
  
!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinput

!----------------------------------------------------------------------
! Fill in the Hamiltonian parameter array
!----------------------------------------------------------------------
  call hampar

!----------------------------------------------------------------------
! Compute the Hamiltonian matrix
!----------------------------------------------------------------------
  call calc_hamiltonian

!----------------------------------------------------------------------
! Diagonalise the Hamiltonian matrix
!----------------------------------------------------------------------
  call diag_hamiltonian
  
contains

!######################################################################

  subroutine initialise

    use global
    
    implicit none

!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    ! No. grid points
    Npts=0

    ! Grid bounds
    bounds=-999.0d0

    ! Name of the model Hamiltonian
    aham=''
    
    return
    
  end subroutine initialise
    
!######################################################################

  subroutine rdinput

    use global
    
    implicit none

    integer           :: i
    character(len=80) :: string1,string2

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    i=0

5   i=i+1
    
    call getarg(i,string1)

    if (string1.eq.'-grid') then
       ! Lower bound of the grid
       i=i+1
       call getarg(i,string2)
       read(string2,*) bounds(1)
       ! Upper bound of the grid
       i=i+1
       call getarg(i,string2)
       read(string2,*) bounds(2)
       ! No. grid points
       i=i+1
       call getarg(i,string2)
       read(string2,*) Npts
       ! Make sure that the no. grid points is odd
       if (mod(Npts,2).eq.0) Npts=Npts+1
       ! Grid spacing
       dx=(bounds(2)-bounds(1))/(Npts-1)       
    else if (string1.eq.'-ham') then
       ! Hamiltonian name
       i=i+1
       call getarg(i,aham)
       ! Set the Hamiltonian number
       if (aham.eq.'morsebk') then
          iham=1
       else if (aham.eq.'morseg') then
          iham=2
       else
          write(6,'(/,2x,a,/)') 'Unknown Hamiltonian: '//trim(aham)
       stop
       endif
    else
       write(6,'(/,2x,a,/)') 'Unknown keyword: '//trim(string1)
       stop
    endif

    if (i.lt.iargc()) goto 5

!----------------------------------------------------------------------
! Check the input
!----------------------------------------------------------------------
    ! Grid information
    if (Npts.eq.0) then
       write(6,'(/,2x,a,/)') 'The grid information has not been given'
       stop
    endif

    ! Name of the model potential
    if (aham.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the model Hamiltonian has &
            not been given'
       stop
    endif

    return
    
  end subroutine rdinput

!######################################################################

  subroutine hampar

    use global
    
    implicit none

!----------------------------------------------------------------------
! Number of parameters
!----------------------------------------------------------------------
    select case(iham)

    case(1) ! Balint-Kurti Morse potential
       nvpar=3

    case(2) ! Morse potential used by Guo in
            ! J. Chem. Phys., 105, 1311 (1996)
       nvpar=3

    end select
       
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(vpar(nvpar))
    vpar=0.0d0

!----------------------------------------------------------------------
! Fill in the Hamiltonian parameter array
!----------------------------------------------------------------------
    select case(iham)

    case(1) ! Balint-Kurti Morse potential
       ! D
       vpar(1)=0.1744d0
       ! beta
       vpar(2)=1.02764d0
       ! x_e
       vpar(3)=1.40201d0
       ! reduced mass of H2
       mass=1822.888486192d0/2.0d0

    case(2) ! Guo Morse potential
       ! Prefactor
       vpar(1)=1000.0d0
       ! Exponential 1
       vpar(2)=0.3d0
       ! Exponential 2
       vpar(3)=0.15d0
       ! Mass
       mass=1.0d0

    end select
       
    return
    
  end subroutine hampar

!######################################################################

  subroutine calc_hamiltonian

    use global
    
    implicit none

    integer                         :: i,j,l
    integer                         :: n
    real(dp)                        :: xi,pi,ftmp
    real(dp), dimension((Npts-1)/2) :: Tl

    pi=2.0d0*acos(0.0d0)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hmat(Npts,Npts))
    hmat=0.0d0
    
!----------------------------------------------------------------------
! Potential contribution
!----------------------------------------------------------------------
    do i=1,Npts
       xi=bounds(1)+(i-1)*dx
       hmat(i,i)=potfunc(xi)
    enddo
    
!----------------------------------------------------------------------
! Kinetic energy operator contribution
!----------------------------------------------------------------------
    n=(Npts-1)/2
    
    ! Precalculation of the T_l terms
    do l=1,n
       Tl(l)=(2.0d0/mass)*(pi*l/(Npts*dx))**2
    enddo

    ! Kinetic energy contribution
    do i=1,Npts
       do j=i,Npts

          ftmp=0.0d0
          do l=1,n
             ftmp=ftmp+cos(l*2.0d0*pi*(i-j)/Npts)*Tl(l)
          enddo
          hmat(i,j)=hmat(i,j)+ftmp*2.0d0/Npts

          hmat(j,i)=hmat(i,j)
       enddo
    enddo
          
    return
    
  end subroutine calc_hamiltonian
    
!######################################################################

  function potfunc(x)

    use global
    
    implicit none

    real(dp) :: x,potfunc,V

    select case(iham)
       
    case(1) ! Balint-Kurti Morse potential
       call bkmorse(x,V)

    case(2) ! Guo Morse potential
       call guomorse(x,V)

    end select

    potfunc=V
    
    return
    
  end function potfunc

!######################################################################

  subroutine bkmorse(x,V)

    use global
    
    implicit none

    real(dp) :: x,V,D,beta,xe

    D=vpar(1)
    beta=vpar(2)
    xe=vpar(3)

    V=D*(1.0d0-exp(-beta*(x-xe)))**2
    
    return
    
  end subroutine bkmorse

!######################################################################

  subroutine guomorse(x,V)

    use global

    implicit none

    real(dp) :: x,V,p1,p2,p3

    p1=vpar(1)
    p2=vpar(2)
    p3=vpar(3)
    
    V=p1*(exp(-p2*x)-2.0d0*exp(-p3*x))

    if (V.gt.100.0d0) V=100.0d0
    
    return
    
  end subroutine guomorse
    
!######################################################################

  subroutine diag_hamiltonian

    use global
    
    implicit none

    integer               :: workdim,error,i
    real(dp), allocatable :: work(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    workdim=3*Npts
    allocate(work(workdim))

    allocate(eigvec(Npts,Npts))
    eigvec=0.0d0

    allocate(eigval(Npts))
    eigval=0.0d0
    
!----------------------------------------------------------------------
! Diagonalise the Hamiltonian matrix
!----------------------------------------------------------------------
    eigvec=hmat
    
    call dsyev('V','U',Npts,eigvec,Npts,eigval,work,workdim,error)

    if (error.ne.0) then
       write(6,'(/,2x,a)') 'Error in the diagonalisation of H'
       stop
    endif

    do i=1,Npts
       write(6,'(i5,2x,F11.6)') (i-1),eigval(i)
    enddo

!    do i=1,Npts
!       print*,(i-1)*dx,eigvec(i,16)
!    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(work)
    
    return
    
  end subroutine diag_hamiltonian
    
!######################################################################
  
end program fgh1d
