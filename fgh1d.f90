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
  real(dp), dimension(2) :: gbounds

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

  ! Chebyshev order-domain autocorrelation function
  integer                :: chebyord
  real(dp), dimension(2) :: sbounds
  real(dp), allocatable  :: auto(:)
  real(dp), allocatable  :: hnorm(:,:)
  logical                :: lcheby
  
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

!----------------------------------------------------------------------
! Output the eigenvalues
!----------------------------------------------------------------------
  call wreig
  
!----------------------------------------------------------------------
! Calculation of the Chebyshev order-domain autocorrelation function
!----------------------------------------------------------------------
  if (lcheby) then
     ! Calculate the order-domain autocorrelation function
     call chebyauto
     ! Output the order-domain autocorrelation function
     call wrauto     
  endif
     
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
    gbounds=-999.0d0

    ! Name of the model Hamiltonian
    aham=''

    ! Chebyshev order-domain autocorrelation function calculation
    lcheby=.false.
    chebyord=0    
    
    return
    
  end subroutine initialise
    
!######################################################################

  subroutine rdinput

    use global
    
    implicit none

    integer           :: i
    character(len=80) :: string1,string2

!----------------------------------------------------------------------
! Exit if no arguments were given
!----------------------------------------------------------------------
    if (iargc().eq.0) then
       write(6,'(/,2x,a,/)') 'No arguments given'
       stop
    endif

!----------------------------------------------------------------------
! Print the input keywords if the -h flag was given
!----------------------------------------------------------------------
    call getarg(1,string1)
    if (string1.eq.'-h') call prhelp
    
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
       read(string2,*) gbounds(1)
       ! Upper bound of the grid
       i=i+1
       call getarg(i,string2)
       read(string2,*) gbounds(2)
       ! No. grid points
       i=i+1
       call getarg(i,string2)
       read(string2,*) Npts
       ! Make sure that the no. grid points is odd
       if (mod(Npts,2).eq.0) Npts=Npts+1
       ! Grid spacing
       dx=(gbounds(2)-gbounds(1))/(Npts-1)       
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
    else if (string1.eq.'-cheby') then
       lcheby=.true.
       i=i+1
       call getarg(i,string2)
       read(string2,*) chebyord
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

  subroutine prhelp

    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Write the input options to screen then exit
!----------------------------------------------------------------------
    ! Purpose
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)') 'Purpose'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)') 'Calculates the eigenpairs of a 1D model &
         potential using the Fourier Grid Hamiltonian method.'
    write(6,'(a)') 'Additionally, the Chebyshev order-domain &
         autocorrelation function can be calculated.'
    
    ! Usage
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)') 'Usage'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)') 'fgh1d -grid -ham (-cheby)'

    ! Options
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)')   'Options'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)')     '-ham hamname            : &
         The name of the model Hamiltonian is ''hamname''.'
    write(6,'(a)')     '-grid Ea Eb N           : &
         The grid has N points evenly distributed in the interval &
         [Ea,Eb]'
     write(6,'(a)')     '-cheby K                : &
         The Chebyshev order-domain autocorrelation function &
         will be calculated to order K'

     stop
     
    return
    
  end subroutine prhelp
    
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
       xi=gbounds(1)+(i-1)*dx
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
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(work)
    
    return
    
  end subroutine diag_hamiltonian

!######################################################################

  subroutine wreig

    use global

    implicit none

    integer           :: unit,i
    character(len=80) :: filename
    
!----------------------------------------------------------------------
! Open the energy file
!----------------------------------------------------------------------
    unit=20
    filename=trim(aham)//'_ener.dat'
    open(unit,file=filename,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the energies to file
!----------------------------------------------------------------------
    do i=1,Npts
       write(unit,'(x,i5,2x,F13.8)') (i-1),eigval(i)
    enddo

!----------------------------------------------------------------------
! Close the energy file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine wreig
    
!######################################################################

  subroutine chebyauto

    use global

    implicit none

    integer               :: i,k
    real(dp), allocatable :: q0(:),qk(:),qkm1(:),qkm2(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(auto(0:2*chebyord))
    auto=0.0d0

    allocate(q0(Npts))
    q0=0.0d0

    allocate(qk(Npts))
    qk=0.0d0

    allocate(qkm1(Npts))
    qkm1=0.0d0

    allocate(qkm2(Npts))
    qkm2=0.0d0
    
!----------------------------------------------------------------------
! Set the spectral bounds
!----------------------------------------------------------------------
    sbounds(1)=1.01d0*eigval(1)
    sbounds(2)=1.01d0*eigval(Npts)

!----------------------------------------------------------------------
! Set up initial vector
!----------------------------------------------------------------------
!    do i=1,Npts
!       call random_number(q0(i))
!    enddo
!    q0=q0/sqrt(dot_product(q0,q0))

    ! TEST
    q0(:)=sum(eigvec(:,1:Npts))
    q0=q0/sqrt(dot_product(q0,q0))
    ! TEST
    
!----------------------------------------------------------------------
! Comput the normalised Hamiltonian matrix
!----------------------------------------------------------------------
    call gethnorm
    
!----------------------------------------------------------------------
! C_0
!----------------------------------------------------------------------
    auto(0)=dot_product(q0,q0)

!----------------------------------------------------------------------
! Calculate the Chebyshev order-domain autocorrelation function
!----------------------------------------------------------------------
    ! Initialisation
    qkm1=q0

    ! Loop over Chebyshev polynomials of order k >= 1
    do k=1,chebyord

       ! Calculate the kth Chebyhev polynomial-vector product
       !
       !  H_norm * q_k-1
       qk=matmul(hnorm,qkm1)
       !
       ! 2 * H_norm * q_k-1
       if (k.gt.1) qk=2.0d0*qk
       !
       ! q_k = 2 * H_norm * q_k-1 - q_k-2 (k>2)
       if (k.gt.1) qk=qk-qkm2
       
       ! Calculate C_k
       auto(k)=dot_product(q0,qk)

       ! Calculate C_2k and C_2k-1
       if (k.gt.chebyord/2) then
          auto(2*k)=2.0d0*dot_product(qk,qk)-auto(0)
          auto(2*k-1)=2.0d0*dot_product(qkm1,qk)-auto(1)
       endif

       ! Update qkm1 and qkm2
       qkm2=qkm1       
       qkm1=qk
       qk=0.0d0
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(q0)
    deallocate(qk)
    deallocate(qkm1)
    deallocate(qkm2)
    
    return
    
  end subroutine chebyauto

!######################################################################

  subroutine gethnorm

    use global

    implicit none

    integer  :: i
    real(dp) :: DeltaE,Emin

!    ! TEST
!    integer :: error
!    real(dp), dimension(3*Npts)    :: work
!    real(dp), dimension(Npts,Npts) :: vec
!    real(dp), dimension(Npts)      :: val
!    ! TEST
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(hnorm(Npts,Npts))
    hnorm=0.0d0
    
!----------------------------------------------------------------------
! Compute the normalised Hamiltonian matrix
!----------------------------------------------------------------------
    Emin=sbounds(1)
    DeltaE=sbounds(2)-sbounds(1)

!    ! TEST
!    Emin=eigval(1)
!    DeltaE=eigval(Npts)-eigval(1)
!    ! TEST
    
    hnorm=hmat
    do i=1,Npts
       hnorm(i,i)=hnorm(i,i)-0.5d0*DeltaE-Emin
    enddo
    hnorm=2.0d0*hnorm/DeltaE


!    ! TEST
!    vec=hnorm    
!    call dsyev('V','U',Npts,vec,Npts,val,work,3*Npts,error)
!    if (error.ne.0) then
!       write(6,'(/,2x,a)') 'Error in the diagonalisation of Hnorm'
!       stop
!    endif
!    do i=1,Npts
!       print*,i,val(i)
!    enddo
!    stop
!    ! TEST

    
    return
    
  end subroutine gethnorm

!######################################################################

  subroutine wrauto

    use global

    implicit none

    integer :: unit,k
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    unit=30
    open(unit,file='chebyauto',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the file header
!----------------------------------------------------------------------
    write(unit,'(a,2(2x,E21.14),/)') '#    Spectral bounds:',&
         sbounds(1),sbounds(2)
    write(unit,'(a)')   '#    Order [k]    C_k'
    
!----------------------------------------------------------------------
! Write the order-domain autocorrelation function to file
!----------------------------------------------------------------------
    do k=0,chebyord*2
       write(unit,'(i6,11x,E21.14)') k,auto(k)
    enddo
    
!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine wrauto
    
!######################################################################
  
end program fgh1d
