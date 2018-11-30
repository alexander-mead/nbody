PROGRAM nbody

  USE constants
  USE string_operations
  USE file_info
  USE vectors
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: xi(:,:), vi(:,:), m(:), x(:,:), v(:,:), xnew(:,:), vnew(:,:), xres(:,:,:), tres(:), vres(:,:,:)
  INTEGER :: i, j, k, np, it
  INTEGER :: icm, idim  
  REAL :: tf, dt, e_init, e_final, xcm(3), vcm(3), l_init(3), l_final(3), e, enew, de, t, dtmax, eratio
  REAL :: L(3), Lnew(3), Lmod, Lmodnew, Lratio, dL, acc
  CHARACTER(len=256) :: infile, time, boost, accuracy, directory, dimens
  LOGICAL :: mr_logic, iE, iL
  
  INTEGER, PARAMETER :: iforce=1    ! Change the potential function
  REAL, PARAMETER :: G=1.           ! Gravitational constant
  REAL, PARAMETER :: ti=0.          ! Initial time (no reason to not be zero)
  REAL, PARAMETER :: soft=0.        ! Gravitational softening [AU; need to set iforce=7]
  INTEGER, PARAMETER :: n=1001      ! Number of time-steps to output (linear in t; need to have 101, 1001 etc.)
  REAL, PARAMETER :: min_cons=1e-12 ! If the conserved quanities are below this number then their conservation is not required

  CALL get_command_argument(1,infile)
  IF(infile=='') STOP 'You need to specify an input file'

  CALL get_command_argument(2,time)
  IF(time=='') STOP 'You need to specify a simulation time length (yrs)'
  READ(time,*) tf

  CALL get_command_argument(3,dimens)
  IF(dimens=='') STOP 'You need to specify the number of dimensions (2 or 3)'
  READ(dimens,*) idim
  IF((idim .NE. 2) .AND. (idim .NE. 3)) STOP 'You need to specify the number of dimensions (2 or 3)'

  CALL get_command_argument(4,boost)
  IF(boost=='') STOP 'You need to specify CM boost (1) or not (0)'
  READ(boost,*) icm
  IF((icm .NE. 1) .AND. (icm .NE. 0)) STOP 'You need to specify CM boost (1) or not (0)'

  CALL get_command_argument(5,accuracy)
  IF(accuracy=='') STOP 'You need to specify an accuracy parameter (e.g. 1e-6)'
  READ(accuracy,*) acc

  CALL get_command_argument(6,directory)
  IF(directory=='') STOP 'You need to specify an output directory (e.g. ''data'')'

  INQUIRE(file=infile,exist=mr_logic)
  IF(mr_logic .EQV. .FALSE.) STOP 'Input file does not exist'
  
  WRITE(*,*)
  WRITE(*,*) 'N-body integration code'
  WRITE(*,*) '======================='
  WRITE(*,*)

  ! Read in the initial mass, position and velocity of each particle to be simulated
  CALL read_input(m,xi,vi,np,infile)

  ! I included this because I thought that 2D simulations would be much quicker than 3D
  ! This seems not to be the case though, and the speed-up is very, very marginal
  ! I'm not quite sure why though...
  WRITE(*,*) 'Number of dimensions:', idim
  WRITE(*,*)

  IF(icm==1) THEN
     WRITE(*,*) 'Boosting to CM position'
     xcm=centre_of_mass(np,m,xi)
     vcm=centre_of_mass(np,m,vi)
     DO i=1,np
        DO j=1,idim
           xi(j,i)=xi(j,i)-xcm(j)
           vi(j,i)=vi(j,i)-vcm(j)
        END DO
     END DO
     WRITE(*,*) 'CM move complete'
     WRITE(*,*)
  END IF

  ! Write useful information to the screen
  WRITE(*,*) 'Mass units [Msun]'
  WRITE(*,*) 'Length units [AU]'
  WRITE(*,*) 'Velocity units [2*pi*AU/year]'
  Write(*,*) 'Start time of [years]:', ti
  WRITE(*,*) 'End time of [years]:', tf
  WRITE(*,*) 'Duration of simulation [years]:', tf-ti
  WRITE(*,*) 'Gravitational softening [AU]:', soft
  WRITE(*,*)

  ! Convert to internal time units (Done so that G=1 in the code, rather than 4*pi^2)
  tf=tf*2.*pi

  ALLOCATE(x(3,np),xnew(3,np),v(3,np),vnew(3,np),xres(3,np,n),vres(3,np,n),tres(n))

  ! Set variables to zero
  x=0.
  xnew=0.
  v=0.
  vnew=0.
  t=ti

  ! Write the ICs to be the first line in the results file
  DO k=1,3
     DO i=1,np
        x(k,i)=xi(k,i)
        v(k,i)=vi(k,i)
        xres(k,i,1)=x(k,i)
        vres(k,i,1)=v(k,i)
     END DO
  END DO

  ! Assume checking AM and energy conservation
  iL=.TRUE.
  iE=.TRUE.

  L=angular_momentum(m,x,v)
  Lmod=modulus(L,3)
  e=total_energy(m,xi,vi)

  WRITE(*,*) 'Initial energy:', e
  WRITE(*,*) 'Initial angular momentum:', Lmod
  WRITE(*,*)
  
  ! Don't do AM conservation if |L|=0. or energy conservation if this is zero
  IF(ABS(e)<=min_cons) iE=.FALSE.
  IF(Lmod<=min_cons)   iL=.FALSE.

  ! Set the physical time-step length based on the desired number (n) outputs
  ! Actually there will be n-1 further outputs because n=1 is taken up with ICs
  dt=(tf-ti)/float(n-1)

  ! Maximum timestep is then set
  dtmax=dt

  ! Time stepping integer; dt=dtmax/(2**it)
  it=0

  WRITE(*,*) 'Starting intergation'
  IF(iE) WRITE(*,*) 'Integrating while checking energy conservation'
  IF(iL) WRITE(*,*) 'Integrating while checking angular momentum conservation'
  WRITE(*,*) 'Accuracy parameter:', acc ! User set and is the degree to which AM and energy conservation occur
  WRITE(*,*)

  i=1 ! Must set i=1 here to make calculation work properly
  DO

     ! Calculate energy and AM at the beginning of the step
     IF(iE) THEN
        e=total_energy(m,x,v)
     END IF
     IF(iL) THEN
        L=angular_momentum(m,x,v)
        Lmod=modulus(L,3)
     END IF

     ! This integrates the system from t to t+dt
     CALL rk4(np,x(:,:),xnew(:,:),v(:,:),vnew(:,:),dt)

     ! Calculate energy and AM at the end of the step
     IF(iE) THEN
        enew=total_energy(m,xnew,vnew)
     END IF
     IF(iL) THEN
        Lnew=angular_momentum(m,x,v)
        Lmodnew=modulus(Lnew,3)
     END IF

     ! Calculate ratio of before and after energy
     IF(iE) THEN
        eratio=ABS(-1+enew/e)
        de=acc*dt/tf
     END IF

     ! Calculate ratio of before and after angular momentum
     IF(iL) THEN
        Lratio=ABS(-1+Lmodnew/Lmod)
        dL=acc*dt/tf
     END IF

     IF(((iE .EQV. .FALSE.) .OR. (eratio<de)) .AND. ((iL .EQV. .FALSE.) .OR. (Lratio<dL))) THEN
        
        ! In this case update the positions and move on to the next step
        x=xnew
        v=vnew
        !t=t+dt ! NB. I think errors accrew in the time calculation (big sum) this was so that tf won't be exactly what you specify

        ! Although this might not be the case because each timestep is a power of two thing 1/2**N
        t=t+dt

        ! This is only approximately accurate for output time, but is easiest
        ! As set there is no guarenee that the output time will be an integer multiple of dtmax because this exact value
        ! can be missed. I can't think of an easy way to fix this at the moment.
        ! Still, the *value* of time should be accurate
        IF(t>=dtmax*float(i)) THEN

           ! Using these integer routines should be 100% accurate for output time
           ! but seems to be super slow, must be all the integer operations
           ! Actually this will still miss doing integer multiples of 'dt_{max}' (e.g. 5/8+1/2)
           !yint(2)=it
           !xint=add2powfrac(xint,yint)
           !IF(xint(2)==0) THEN
           !xint(1)=0
           !xint(2)=1
           
           ! Tells the user how long has elapsed in terms of time
           IF(i==CEILING(0.1*DBLE(n))) WRITE(*,*) '10% Complete'
           IF(i==CEILING(0.2*DBLE(n))) WRITE(*,*) '20% Complete'
           IF(i==CEILING(0.3*DBLE(n))) WRITE(*,*) '30% Complete'
           IF(i==CEILING(0.4*DBLE(n))) WRITE(*,*) '40% Complete'
           IF(i==CEILING(0.5*DBLE(n))) WRITE(*,*) '50% Complete'
           IF(i==CEILING(0.6*DBLE(n))) WRITE(*,*) '60% Complete'
           IF(i==CEILING(0.7*DBLE(n))) WRITE(*,*) '70% Complete'
           IF(i==CEILING(0.8*DBLE(n))) WRITE(*,*) '80% Complete'
           IF(i==CEILING(0.9*DBLE(n))) WRITE(*,*) '90% Complete'
           i=i+1 ! i will be 1+1=2 on first pass so this does not overwrite IC
           xres(:,:,i)=x(:,:)
           vres(:,:,i)=v(:,:)
           tres(i)=t
           IF(i==n) THEN
              WRITE(*,*) '100% Complete'
              EXIT
           END IF
        END IF
        IF(it .NE. 0) THEN
           ! Try an increase the time-step as long as it below the maximum value
           !dt=dt*2
           it=it-1
           dt=dtmax/float(2**it)
        END IF
        
     ELSE

        ! Otherwise decrease the time-step and do the calculation again
        it=it+1
        dt=dtmax/float(2**it)

        ! This is an attempt to make the code exit if the timestep becomes too small
        ! I've not really checked to see what an 'optimum' value is here
        IF(it==30) THEN
           WRITE(*,*) 'Error - collision detected exiting'
           WRITE(*,*) 'Collision at time:', t/(2.*pi)
           WRITE(*,*) 'Finishing time:', tf/(2.*pi)
           DO j=i+1,n
              xres(:,:,j)=x(:,:)
              vres(:,:,j)=v(:,:)
              tres(j)=0
           END DO
           EXIT
        END IF

     END IF

  END DO
  WRITE(*,*)
  ! Integration is now finished
  
  ! Calculate initial and final energies, AM and CM and VCM
  e_init=total_energy(m,xi(:,:),vi(:,:))
  e_final=total_energy(m,x(:,:),v(:,:))
  l_init=angular_momentum(m,xi(:,:),vi(:,:))
  l_final=angular_momentum(m,x(:,:),v(:,:))
  xcm=centre_of_mass(np,m,x(:,:))
  vcm=centre_of_mass(np,m,v(:,:))

  WRITE(*,*) 'Final CM position [AU]'
  WRITE(*,*) xcm(1)
  WRITE(*,*) xcm(2)
  WRITE(*,*) xcm(3)
  WRITE(*,*)

  WRITE(*,*) 'Final CM velocity [2*pi*AU/year]'
  WRITE(*,*) vcm(1)
  WRITE(*,*) vcm(2)
  WRITE(*,*) vcm(3)
  WRITE(*,*)

  ! Angular momentum in x
  IF(l_init(1) .NE. 0.) THEN
     WRITE(*,*) 'Initial Lx:', l_init(1)
     WRITE(*,*) 'Final Lx:', l_final(1)
     WRITE(*,*) 'Ratio:', l_final(1)/l_init(1)
     WRITE(*,*)
  END IF

  ! Angular momentum in y
  IF(l_init(2) .NE. 0.) THEN
     WRITE(*,*) 'Initial Ly:', l_init(2)
     WRITE(*,*) 'Final Ly:', l_final(2)
     WRITE(*,*) 'Ratio:', l_final(2)/l_init(2)
     WRITE(*,*)
  END IF

  ! Angular momentum in z
  IF(l_init(3) .NE. 0.) THEN
     WRITE(*,*) 'Initial Lz:', l_init(3)
     WRITE(*,*) 'Final Lz:', l_final(3)
     WRITE(*,*) 'Ratio:', l_final(3)/l_init(3)
     WRITE(*,*)
  END IF

  ! Total angular momentum
  IF(e_init .NE. 0.) THEN
     WRITE(*,*) 'Initial energy:', e_init
     WRITE(*,*) 'Final energy:', e_final
     WRITE(*,*) 'Ratio:', e_final/e_init
     WRITE(*,*)
  END IF

  ! This the writes out the results, it does this only when the calculation is complete
  CALL write_results(np,n,tres,xres,vres,m,directory)

CONTAINS

  FUNCTION centre_of_mass(n,m,x)

    ! Calculate the CM vector from all particles
    IMPLICIT NONE
    REAL, INTENT(IN) :: m(:), x(:,:)
    REAL :: centre_of_mass(3)
    INTEGER :: i, j, n

    ! Set sum variable to zero
    centre_of_mass=0.

    ! Sum over dimensions and numbers of particles
    DO j=1,3
       DO i=1,n
          centre_of_mass(j)=centre_of_mass(j)+m(i)*x(j,i)
       END DO
    END DO

    centre_of_mass=centre_of_mass/SUM(m)

  END FUNCTION centre_of_mass

  FUNCTION total_energy(m,x,v)

    ! Calculates the total energy of all particles
    IMPLICIT NONE
    INTEGER :: n, i, j
    REAL :: kin, pot, total_energy, d
    REAL, INTENT(IN) :: m(:), x(:,:), v(:,:)

    ! Set the sum variables to zero
    kin=0.
    pot=0.

    n=SIZE(x,2)

    ! Calculate kinetic energy
    DO i=1,n
       kin=kin+m(i)*((v(1,i)**2.)+(v(2,i)**2.)+(v(3,i)**2.))/2.
    END DO

    ! Calculate potential energy
    DO i=1,n
       DO j=1,i-1
          d=distance(x(:,i),x(:,j),idim)
          pot=pot+G*m(i)*m(j)*gravitational_potential(d)
       END DO
    END DO

    ! Total is sum
    total_energy=kin+pot

  END FUNCTION total_energy

  FUNCTION angular_momentum(m,x,v)

    ! This calculates the total angular momentum (about the origin) of the particles L=sum r x p
    IMPLICIT NONE
    REAL :: angular_momentum(3)
    REAL, INTENT(IN) :: m(:), x(:,:), v(:,:)
    INTEGER :: i, n

    ! Set sum variable to zero
    angular_momentum=0.

    n=SIZE(x,2)

    DO i=1,n
       angular_momentum=angular_momentum+m(i)*cross_product(x(:,i),v(:,i))
    END DO

  END FUNCTION angular_momentum

  FUNCTION gravitational_potential(r)

    ! Potential energy without the G*m1*m2 pre-factor
    IMPLICIT NONE
    REAL :: gravitational_potential
    REAL, INTENT(IN) :: r
    REAL :: rr    

    ! Softening length
    rr=r+soft

    ! Different central potentials
    IF(iforce==1) THEN
       gravitational_potential=-1./rr
    ELSE IF(iforce==2) THEN
       gravitational_potential=-1./rr**2
    ELSE IF(iforce==3) THEN
       gravitational_potential=rr+rr**(-1)
    ELSE IF(iforce==4) THEN
       gravitational_potential=rr**2.+rr**(-2)
    ELSE IF(iforce==5) THEN
       gravitational_potential=rr**0.5+rr**(-0.5)
    ELSE IF(iforce==6) THEN
       gravitational_potential=rr**(-3)-rr
    ELSE
       STOP 'GRAVITATIONAL_POTENTIAL: Error, iforce not specified correctly'
    END IF
    
  END FUNCTION gravitational_potential

  FUNCTION gravitational_force(r)

    ! Calculates the central force without the G*m1*m2 pre-factor
    IMPLICIT NONE
    REAL :: gravitational_force
    REAL, INTENT(IN) :: r
    REAL :: rr

    ! Softening length
    rr=r+soft

    ! Force = -Grad V(r) (NB. *MINUS* grad V); V is a function of r only
    IF(iforce==1) THEN
       gravitational_force=-1./rr**2
    ELSE IF(iforce==2) THEN
       gravitational_force=-2./rr**3
    ELSE IF(iforce==3) THEN
       gravitational_force=-1.+rr**(-2)
    ELSE IF(iforce==4) THEN
       gravitational_force=-2.*rr+2*rr**(-3)
    ELSE IF(iforce==5) THEN
       gravitational_force=-0.5*rr**(-0.5)+0.5*rr**(-1.5)
    ELSE IF(iforce==6) THEN
       gravitational_force=3.*rr**(-4)-rr**(-2)
    ELSE
       STOP 'GRAVITATIONAL_FORCE: Error, iforce not specified correctly'
    END IF
    
  END FUNCTION gravitational_force

  SUBROUTINE force_matrix(n,x,m,F)

    ! This is the force matrix calculation does the i .NE. j forces and reflects along the diagonal
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(:,:), m(:)
    INTEGER :: n
    REAL, INTENT(OUT) :: F(:,:,:)
    INTEGER :: i, j, k
    REAL :: d   

    ! The diagonal is never investigated so will remain zero (anti-symmetry, no self force)
    F=0.

    DO j=1,n
       DO i=1,j-1

          ! Calculate the particle-particle distances depending on the number of dimensions
          d=distance(x(:,i),x(:,j),idim)

          ! Compute all elements of the force matrix, anti-symmetry is enforced; radial vector comes in at the end
          DO k=1,idim
             F(k,i,j)=G*m(i)*m(j)*gravitational_force(d)*(x(k,i)-x(k,j))/d
             F(k,j,i)=-F(k,i,j)
          END DO
          
       END DO
    END DO

  END SUBROUTINE force_matrix

  SUBROUTINE rk4(n,xin,xout,vin,vout,dt)

    IMPLICIT NONE
    REAL, INTENT(IN) :: xin(:,:), vin(:,:), dt
    REAL, INTENT(OUT) :: xout(:,:), vout(:,:)
    INTEGER, INTENT(IN) :: n
    REAL, ALLOCATABLE :: x(:,:), v(:,:), F(:,:,:)!, disp(:,:,:)
    REAL, ALLOCATABLE :: kx1(:,:), kx2(:,:), kx3(:,:), kx4(:,:)
    REAL, ALLOCATABLE :: kv1(:,:), kv2(:,:), kv3(:,:), kv4(:,:)
    REAL :: accl
    INTEGER :: i, j, k

    ! This is a generic routine to carry out the RK4 algorithm between points t and t+dt'
    ! for n bodies experiencing forces between each other

    ! Allocate arrays
    ALLOCATE(kx1(idim,n),kx2(idim,n),kx3(idim,n),kx4(idim,n))
    ALLOCATE(kv1(idim,n),kv2(idim,n),kv3(idim,n),kv4(idim,n))
    ALLOCATE(x(idim,n),v(idim,n),F(idim,n,n))

    ! This means that x, v in this subroutine can be either 2 or 3 dimensions
    DO k=1,idim
       DO i=1,n
          x(k,i)=xin(k,i)
          v(k,i)=vin(k,i)
       END DO
    END DO
    kv1=0.
    kv2=0.
    kv3=0.
    kv4=0.

    kx1=0.
    kx2=0.
    kx3=0.
    kx4=0.

    ! Step 1
    ! Calculates the force matrix
    CALL force_matrix(n,x,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                !accl=G*m(j)*disp(k,i,j)
                accl=F(k,i,j)/m(i)
                kv1(k,i)=kv1(k,i)+accl*dt
             END IF
          END DO
          kx1(k,i)=v(k,i)*dt
       END DO
    END DO

    ! Step 2
    CALL force_matrix(n,x+kx1/2.,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=F(k,i,j)/m(i)
                kv2(k,i)=kv2(k,i)+accl*dt
             END IF
          END DO
          kx2(k,i)=(v(k,i)+kv1(k,i)/2.)*dt
       END DO
    END DO

    ! Step 3
    CALL force_matrix(n,x+kx2/2.,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=F(k,i,j)/m(i)
                kv3(k,i)=kv3(k,i)+accl*dt
             END IF
          END DO
          kx3(k,i)=(v(k,i)+kv2(k,i)/2.)*dt
       END DO
    END DO

    ! Step 4
    CALL force_matrix(n,x+kx3,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=F(k,i,j)/m(i)
                kv4(k,i)=kv4(k,i)+accl*dt
             END IF
          END DO
          kx4(k,i)=(v(k,i)+kv3(k,i))*dt
       END DO
    END DO

    DEALLOCATE(x,v,F)

    ! Now compute the RK4 result using the standard weighting   
    ! Only run over the number of dimensions
    DO k=1,idim
       DO i=1,n
          xout(k,i)=xin(k,i)+(kx1(k,i)+(2.*kx2(k,i))+(2.*kx3(k,i))+kx4(k,i))/6.
          vout(k,i)=vin(k,i)+(kv1(k,i)+(2.*kv2(k,i))+(2.*kv3(k,i))+kv4(k,i))/6.
       END DO
    END DO

    ! Need to set some values for the z component in 2 dimensional case. These should be 0 from before anyway
    IF(idim==2) THEN
       DO i=1,n
          xout(3,i)=xin(3,i)
          vout(3,i)=vin(3,i)
       END DO
    END IF

    DEALLOCATE(kx1,kx2,kx3,kx4)
    DEALLOCATE(kv1,kv2,kv3,kv4)

  END SUBROUTINE rk4

  SUBROUTINE read_input(m,x,v,n,file_name)

    ! Reads the input file
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: m(:), x(:,:), v(:,:)
    CHARACTER(len=256), INTENT(IN) :: file_name
    INTEGER, INTENT(OUT) :: n  

    n=file_length(file_name,verbose=.TRUE.)

    WRITE(*,*) 'Particles in simulation:', n

    ALLOCATE(x(3,n),v(3,n),m(n))

    WRITE(*,*) 'Reading input'

    OPEN(7,file=file_name)
    DO i=1,n
       
       READ(7,*) m(i), x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i)

       ! If only using 2 dimensions then ensure z compoents are set to zero
       IF(idim==2) THEN
          x(3,i)=0.
          v(3,i)=0.
       END IF
       
    END DO
    CLOSE(7)

    WRITE(*,*) 'Finished reading input'
    WRITE(*,*)

  END SUBROUTINE read_input

  SUBROUTINE write_results(np,n,t,x,v,m,dir)

    ! Writes out the results
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, np
    REAL, ALLOCATABLE :: t(:), x(:,:,:), v(:,:,:), m(:)
    CHARACTER(len=256) :: fname, stem, ext, dir
    INTEGER :: i, j

    ! Convert from internal to external time units
    t=t/(2.*pi)
  
    stem=TRIM(dir)//'/particle_'
    ext='.dat'

    ! Writes out the n positions and velocities
    WRITE(*,*) 'Writing output:'
    DO j=1,np
       fname=number_file(stem,j,ext)
       OPEN(7,file=fname)
       WRITE(*,*) TRIM(fname)
       DO i=1,n
          WRITE(7,fmt='(7F20.10)') t(i), x(1,j,i), x(2,j,i), x(3,j,i), v(1,j,i), v(2,j,i), v(3,j,i)
       END DO
       CLOSE(7)
    END DO

    ! Writes out the final positions and velocities in the format of an input file in case the calculation needs to be resumed
    fname=TRIM(dir)//'/end.dat'
    WRITE(*,*) TRIM(fname)
    OPEN(7,file=fname)
    DO i=1,np
       WRITE(7,fmt='(7F15.7)') m(i), x(1,i,n), x(2,i,n), x(3,i,n), v(1,i,n), v(2,i,n), v(3,i,n)
    END DO
    CLOSE(7)
    WRITE(*,*)

  END SUBROUTINE write_results

END PROGRAM nbody
