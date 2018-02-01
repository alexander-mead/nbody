PROGRAM nbody

  IMPLICIT NONE
  REAL*16, ALLOCATABLE :: xi(:,:), vi(:,:), m(:), x(:,:), v(:,:), xnew(:,:), vnew(:,:), xres(:,:,:), tres(:), vres(:,:,:)
  INTEGER :: i, j, k, np, n, it!, xint(2), yint(2)
  INTEGER :: icm, idim  
  REAL*16 :: tf, dt, e_init, e_final, xcm(3), vcm(3), l_init(3), l_final(3), e, enew, de, t, dtmax, eratio
  REAL*16 :: L(3), Lnew(3), Lmod, Lmodnew, Lratio, dL, acc
  CHARACTER(len=256) :: infile, time, boost, accuracy, directory, dimens
  LOGICAL :: mr_logic, iE, iL
  
  INTEGER, PARAMETER :: iforce=1 !Change the potential function
  REAL*16, PARAMETER :: G=1.d0 !Gravitational constant
  REAL*16, PARAMETER :: pi=3.141592654d0 !pi
  REAL*16, PARAMETER :: ti=0.d0 !Initial time (no reason to not be zero)
  REAL*16, PARAMETER :: soft=0.0d0 !Gravitational softening [AU; need to set iforce=7]

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

  !Number of time-steps to output (linear in t)
  n=1001 !Needs to have the '1' at the end to make the per-cent calculation work
  
  WRITE(*,*)
  WRITE(*,*) 'N-body integration code'
  WRITE(*,*) '======================='
  WRITE(*,*)

  !Reads in the initial mass, position and velocity of each particle to be simulated
  CALL read_input(m,xi,vi,np,infile)

  !I included this because I thought that 2D simulations would be much quicker than 3D
  !This seems not to be the case though, and the speed-up is very, very marginal
  !I'm not quite sure why though...
  WRITE(*,*) 'Number of dimensions:', idim
  WRITE(*,*)

  IF(icm==1) THEN
     WRITE(*,*) 'Boosting to CM position'
     xcm=cm(np,m,xi)
     vcm=cm(np,m,vi)
     DO i=1,np
        DO j=1,idim
           xi(j,i)=xi(j,i)-xcm(j)
           vi(j,i)=vi(j,i)-vcm(j)
        END DO
     END DO
     WRITE(*,*) 'CM move complete'
     WRITE(*,*)
  END IF

  WRITE(*,*) 'Mass units [Msun]'
  WRITE(*,*) 'Length units [AU]'
  WRITE(*,*) 'Velocity in [2*pi*AU/year]:'
  WRITE(*,*) 'Start time of [years]:', tf-ti
  WRITE(*,*) 'End time of [years]:', tf-ti
  WRITE(*,*) 'Duration of simulation [years]:', tf-ti
  WRITE(*,*) 'Gravitational softening [AU]:', soft
  WRITE(*,*)

  !Convert to internal time units (Done so that G=1 in the code, rather than 4*pi^2)
  tf=tf*2.d0*pi

  ALLOCATE(x(3,np),xnew(3,np),v(3,np),vnew(3,np),xres(3,np,n),vres(3,np,n),tres(n))

  !Set variables to zero
  x=0.d0
  xnew=0.d0
  v=0.d0
  vnew=0.d0
  t=ti

  !Write the ICs to be the first line in the results file
  DO k=1,3
     DO i=1,np
        x(k,i)=xi(k,i)
        v(k,i)=vi(k,i)
        xres(k,i,1)=x(k,i)
        vres(k,i,1)=v(k,i)
     END DO
  END DO

  !Assume checking AM and energy conservation
  iL=.TRUE.
  iE=.TRUE.

  L=am(m,x,v)
  Lmod=length(L)
  e=energy(m,xi,vi)

  WRITE(*,*) 'Initial energy:', e
  WRITE(*,*) 'Initial angular momentum:', Lmod
  WRITE(*,*)
  
  !Don't do AM conservation if |L|=0. or energy conservation if this is zero
  IF(ABS(e)<=1d-12) iE=.FALSE.
  IF(Lmod<=1d-12) iL=.FALSE.

  !Set the physical time-step length based on the desired number (n) outputs
  !Actually there will be n-1 further outputs because n=1 is taken up with ICs
  dt=(tf-ti)/float(n-1)

  !Maximum timestep is then set
  dtmax=dt

  !Time stepping integer; dt=dtmax/(2**it)
  it=0
!  xint(1)=0
!  xint(2)=1
!  yint(1)=1
!  yint(2)=1

  WRITE(*,*) 'Starting intergation'
  IF(iE) WRITE(*,*) 'Integrating while checking energy conservation'
  IF(iL) WRITE(*,*) 'Integrating while checking angular momentum conservation'
  WRITE(*,*) 'Accuracy parameter:', acc !User set and is the degree to which AM and energy conservation occur
  WRITE(*,*)

  i=1 !Must set i=1 here to make calculation work properly
  DO

     !Calculate energy and AM at the beginning of the step
     IF(iE) THEN
        e=energy(m,x,v)
     END IF
     IF(iL) THEN
        L=am(m,x,v)
        Lmod=length(L)
     END IF

     !This integrates the system from t to t+dt
     CALL rk4(np,x(:,:),xnew(:,:),v(:,:),vnew(:,:),dt)

     !Calculate energy and AM at the end of the step
     IF(iE) THEN
        enew=energy(m,xnew,vnew)
     END IF
     IF(iL) THEN
        Lnew=am(m,x,v)
        Lmodnew=length(Lnew)
     END IF

     !Calculate ratio of before and after energy
     IF(iE) THEN
        eratio=ABS(-1.d0+enew/e)
        de=acc*dt/tf
     END IF

     !Calculate ratio of before and after angular momentum
     IF(iL) THEN
        Lratio=ABS(-1.d0+Lmodnew/Lmod)
        dL=acc*dt/tf
     END IF

     IF(((iE .EQV. .FALSE.) .OR. (eratio<de)) .AND. ((iL .EQV. .FALSE.) .OR. (Lratio<dL))) THEN
        !In this case update the positions and move on to the next step
        x=xnew
        v=vnew
        !t=t+dt !NB. I think errors accrew in the time calculation (big sum) this was so that tf won't be exactly what you specify
        !Although this might not be the case because each timestep is a power of two thing 1/2**N
        t=t+dt

        !This is only approximately accurate for output time, but is easiest
        !As set there is no guarenee that the output time will be an integer multiple of dtmax because this exact value
        !can be missed. I can't think of an easy way to fix this at the moment.
        !Still, the *value* of time should be accurate
        IF(t>=dtmax*float(i)) THEN

           !Using these integer routines should be 100% accurate for output time
           !but seems to be super slow, must be all the integer operations
           !Actually this will still miss doing integer multiples of 'dt_{max}' (e.g. 5/8+1/2)
           !yint(2)=it
           !xint=add2powfrac(xint,yint)
           !IF(xint(2)==0) THEN
           !xint(1)=0
           !xint(2)=1
           
           !Tells the user how long has elapsed in terms of time
           IF(i==CEILING(0.1*DBLE(n))) WRITE(*,*) '10% Complete'
           IF(i==CEILING(0.2*DBLE(n))) WRITE(*,*) '20% Complete'
           IF(i==CEILING(0.3*DBLE(n))) WRITE(*,*) '30% Complete'
           IF(i==CEILING(0.4*DBLE(n))) WRITE(*,*) '40% Complete'
           IF(i==CEILING(0.5*DBLE(n))) WRITE(*,*) '50% Complete'
           IF(i==CEILING(0.6*DBLE(n))) WRITE(*,*) '60% Complete'
           IF(i==CEILING(0.7*DBLE(n))) WRITE(*,*) '70% Complete'
           IF(i==CEILING(0.8*DBLE(n))) WRITE(*,*) '80% Complete'
           IF(i==CEILING(0.9*DBLE(n))) WRITE(*,*) '90% Complete'
           i=i+1 !i will be 1+1=2 on first pass so this does not overwrite IC
           xres(:,:,i)=x(:,:)
           vres(:,:,i)=v(:,:)
           tres(i)=t
           IF(i==n) THEN
              WRITE(*,*) '100% Complete'
              EXIT
           END IF
        END IF
        !        IF(dt<dtmax/1.1d0) THEN
        IF(it .NE. 0) THEN
           !Try an increase the time-step as long as it below the maximum value
        !          dt=dt*2.d0
           it=it-1
           dt=dtmax/float(2**it)
        END IF
        
     ELSE

        !Otherwise decrease the time-step and do the calculation again
        it=it+1
        dt=dtmax/float(2**it)

        !This is an attempt to make the code exit if the timestep becomes too small
        !I've not really checked to see what an 'optimum' value is here
        IF(it==30) THEN
           WRITE(*,*) 'Error - collision detected exiting'
           WRITE(*,*) 'Collision at time:', t/(2.*pi)
           WRITE(*,*) 'Finishing time:', tf/(2.*pi)
           DO j=i+1,n
              xres(:,:,j)=x(:,:)
              vres(:,:,j)=v(:,:)
              tres(j)=0.d0
           END DO
           EXIT
        END IF

     END IF

  END DO
  WRITE(*,*)
  !Integration is now finished
  
  !Calculate initial and final energies, AM and CM and VCM
  e_init=energy(m,xi(:,:),vi(:,:))
  e_final=energy(m,x(:,:),v(:,:))
  l_init=am(m,xi(:,:),vi(:,:))
  l_final=am(m,x(:,:),v(:,:))
  xcm=cm(np,m,x(:,:))
  vcm=cm(np,m,v(:,:))

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

  IF(l_init(1) .NE. 0.d0) THEN
     WRITE(*,*) 'Initial Lx:', l_init(1)
     WRITE(*,*) 'Final Lx:', l_final(1)
     WRITE(*,*) 'Ratio:', l_final(1)/l_init(1)
     WRITE(*,*)
  END IF
  
  IF(l_init(2) .NE. 0.d0) THEN
     WRITE(*,*) 'Initial Ly:', l_init(2)
     WRITE(*,*) 'Final Ly:', l_final(2)
     WRITE(*,*) 'Ratio:', l_final(2)/l_init(2)
     WRITE(*,*)
  END IF
  
  IF(l_init(3) .NE. 0.d0) THEN
     WRITE(*,*) 'Initial Lz:', l_init(3)
     WRITE(*,*) 'Final Lz:', l_final(3)
     WRITE(*,*) 'Ratio:', l_final(3)/l_init(3)
     WRITE(*,*)
  END IF

  IF(e_init .NE. 0.d0) THEN
     WRITE(*,*) 'Initial energy:', e_init
     WRITE(*,*) 'Final energy:', e_final
     WRITE(*,*) 'Ratio:', e_final/e_init
     WRITE(*,*)
  END IF

  !This the writes out the results, it does this only when the calculation is complete
  CALL results(np,n,tres,xres,vres,m,directory)

CONTAINS

!!$  PURE FUNCTION add2powfrac(x,y)
!!$
!!$    IMPLICIT NONE
!!$    INTEGER*8 :: add2powfrac(2)
!!$    INTEGER*8, INTENT(IN) :: x(2), y(2)
!!$    INTEGER*8 :: a, b, n, m, c, q
!!$
!!$    !Adds two fractions together as long as the denominators
!!$    !Are 2^n with n integer
!!$    !e.g fraction x is x=x(1)/(2**x(2)) !TAKE CARE!
!!$
!!$    a=x(1)
!!$    n=x(2)
!!$    
!!$    b=y(1)
!!$    m=y(2)
!!$    
!!$    IF(n==m) THEN
!!$       q=m
!!$       c=a+b
!!$    ELSE IF(n<m) THEN
!!$       q=m
!!$       c=a*2**(m-n)+b
!!$    ELSE IF(n>m) THEN
!!$       q=n
!!$       c=a+b*(2**(n-m))
!!$    END IF
!!$
!!$    DO
!!$       IF(mod(c,2)==0) THEN
!!$          c=c/2
!!$          q=q-1
!!$       ELSE
!!$          EXIT
!!$       END IF
!!$    END DO
!!$
!!$    add2powfrac(1)=c
!!$    add2powfrac(2)=q
!!$    
!!$  END FUNCTION

  FUNCTION am(m,x,v)

    IMPLICIT NONE
    REAL*16 :: am(3)
    REAL*16, INTENT(IN) :: m(:), x(:,:), v(:,:)
    INTEGER :: i, n

    !This calculates the total angular momentum (about the origin) of the particles L=sum r x p

    am=0.d0

    n=SIZE(x,2)

    DO i=1,n
       am=am+m(i)*cross_product(x(:,i),v(:,i))
    END DO

  END FUNCTION am

  SUBROUTINE force_matrix(n,x,m,F)

    IMPLICIT NONE
    REAL*16, INTENT(IN) :: x(:,:), m(:)
    INTEGER :: n
    REAL*16, INTENT(OUT) :: F(:,:,:)
    INTEGER :: i, j, k
    REAL*16 :: d

    !This is the force matrix calculation does the i .NE. j forces and reflects along the diagonal

    !The diagonal is never investigated so will remain zero (anti-symmetry, no self force)
    F=0.d0

    DO j=1,n
       DO i=1,j-1

          !Calculate the particle-particle distances depending on the number of dimensions
          IF(idim==3) THEN
             d=dist_3D(x(:,i),x(:,j))
          ELSE IF(idim==2) THEN
             d=dist_2D(x(:,i),x(:,j))
          END IF

          !Compute all elements of the force matrix, anti-symmetry is enforced
          !radial vector comes in at the end
          DO k=1,idim
             F(k,i,j)=G*m(i)*m(j)*force(d)*(x(k,i)-x(k,j))/d
             F(k,j,i)=-F(k,i,j)
          END DO
          
       END DO
    END DO

  END SUBROUTINE force_matrix

  SUBROUTINE rk4(n,xin,xout,vin,vout,dt)

    IMPLICIT NONE
    REAL*16, INTENT(IN) :: xin(:,:), vin(:,:), dt
    REAL*16, INTENT(OUT) :: xout(:,:), vout(:,:)
    INTEGER, INTENT(IN) :: n
    REAL*16, ALLOCATABLE :: x(:,:), v(:,:), F(:,:,:)!, disp(:,:,:)
    REAL*16, ALLOCATABLE :: kx1(:,:), kx2(:,:), kx3(:,:), kx4(:,:)
    REAL*16, ALLOCATABLE :: kv1(:,:), kv2(:,:), kv3(:,:), kv4(:,:)
    REAL*16 :: accl
    INTEGER :: i, j, k

    !This is a generic routine to carry out the RK4 algorithm between points t and t+dt'
    !for n bodies experiencing forces between each other

    !Allocate arrays
    ALLOCATE(kx1(idim,n),kx2(idim,n),kx3(idim,n),kx4(idim,n))
    ALLOCATE(kv1(idim,n),kv2(idim,n),kv3(idim,n),kv4(idim,n))
    ALLOCATE(x(idim,n),v(idim,n),F(idim,n,n))

    !This means that x, v in this subroutine can be either 2 or 3 dimensions
    DO k=1,idim
       DO i=1,n
          x(k,i)=xin(k,i)
          v(k,i)=vin(k,i)
       END DO
    END DO
    kv1=0.d0
    kv2=0.d0
    kv3=0.d0
    kv4=0.d0

    kx1=0.d0
    kx2=0.d0
    kx3=0.d0
    kx4=0.d0

    !Step 1
    !Calculates the force matrix
    !CALL force_matrix(n,x,disp)
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

    !Step 2
    CALL force_matrix(n,x+kx1/2.d0,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                !accl=G*m(j)*disp(k,i,j)
                accl=F(k,i,j)/m(i)
                kv2(k,i)=kv2(k,i)+accl*dt
             END IF
          END DO
          kx2(k,i)=(v(k,i)+kv1(k,i)/2.d0)*dt
       END DO
    END DO

    !Step 3
    CALL force_matrix(n,x+kx2/2.d0,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                !accl=G*m(j)*disp(k,i,j)
                accl=F(k,i,j)/m(i)
                kv3(k,i)=kv3(k,i)+accl*dt
             END IF
          END DO
          kx3(k,i)=(v(k,i)+kv2(k,i)/2.d0)*dt
       END DO
    END DO

    !Step 4
    CALL force_matrix(n,x+kx3,m,F)
    DO k=1,idim
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                !accl=G*m(j)*disp(k,i,j)
                accl=F(k,i,j)/m(i)
                kv4(k,i)=kv4(k,i)+accl*dt
             END IF
          END DO
          kx4(k,i)=(v(k,i)+kv3(k,i))*dt
       END DO
    END DO

    DEALLOCATE(x,v,F)

    !Now compute the RK4 result using the standard weighting
    
    !Only run over the number of dimensions
    DO k=1,idim
       DO i=1,n
          xout(k,i)=xin(k,i)+(kx1(k,i)+(2.d0*kx2(k,i))+(2.d0*kx3(k,i))+kx4(k,i))/6.d0
          vout(k,i)=vin(k,i)+(kv1(k,i)+(2.d0*kv2(k,i))+(2.d0*kv3(k,i))+kv4(k,i))/6.d0
       END DO
    END DO

    !Need to set some values for the z component in 2 dimensional case. These should be 0 from before anyway
    IF(idim==2) THEN
       DO i=1,n
          xout(3,i)=xin(3,i)
          vout(3,i)=vin(3,i)
       END DO
    END IF

    DEALLOCATE(kx1,kx2,kx3,kx4)
    DEALLOCATE(kv1,kv2,kv3,kv4)

  END SUBROUTINE rk4

  FUNCTION length(x)

    IMPLICIT NONE
    REAL*16 :: length
    REAL*16, INTENT(IN) :: x(3)

    !Calculates the length of a vector

    length=sqrt(x(1)**2.+x(2)**2.+x(3)**2.)

  END FUNCTION length

!!$  FUNCTION dot_product(x,y)
!!$
!!$    IMPLICIT NONE
!!$    REAL*16 :: dot_product
!!$    REAL*16, INTENT(IN) :: x(3), y(3)
!!$
!!$    !Calculates the dot product of two vectors
!!$
!!$    dot_product=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
!!$
!!$  END FUNCTION dot_product

  FUNCTION cross_product(x,y)

    IMPLICIT NONE
    REAL*16 :: cross_product(3)
    REAL*16, INTENT(IN) :: x(3), y(3)

    !Calculates the cross product of two vectors

    cross_product(1)=x(2)*y(3)-x(3)*y(2)
    cross_product(2)=x(3)*y(1)-x(1)*y(3)
    cross_product(3)=x(1)*y(2)-x(2)*y(1)

  END FUNCTION cross_product

  SUBROUTINE results(np,n,t,x,v,m,dir)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, np
    REAL*16, ALLOCATABLE :: t(:), x(:,:,:), v(:,:,:), m(:)
!    CHARACTER(len=1) :: file_num1, part_num1
!    CHARACTER(len=2) :: file_num2, part_num2
!    CHARACTER(len=3) :: file_num3, part_num3
!    CHARACTER(len=4) :: file_num4, part_num4
    CHARACTER(len=256) :: fname, stem, ext, dir
    INTEGER :: i, j

    !Writes out the results

    !Convert from internal to external time units
    t=t/(2.d0*pi)
  
    stem=TRIM(dir)//'/particle'
    ext='.dat'

    !Writes out the n positions and velocities
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

    !Writes out the final positions and velocities in the format of an input file in case the calculation needs to be resumed
    fname=TRIM(dir)//'/end.dat'
    WRITE(*,*) TRIM(fname)
    OPEN(7,file=fname)
    DO i=1,np
       WRITE(7,fmt='(7F15.7)') m(i), x(1,i,n), x(2,i,n), x(3,i,n), v(1,i,n), v(2,i,n), v(3,i,n)
    END DO
    CLOSE(7)
    WRITE(*,*)

  END SUBROUTINE results

  FUNCTION number_file(fbase,i,fext)

    IMPLICIT NONE
    CHARACTER(len=256) number_file, fbase, fext
    CHARACTER(len=4) num4
    CHARACTER(len=3) num3
    CHARACTER(len=2) num2
    CHARACTER(len=1) num1
    INTEGER :: i

    !A general routine for producing files with a numbered extension

    IF(i<10) THEN
       WRITE(num1,fmt='(I1)') i
       number_file=TRIM(fbase)//'_'//TRIM(num1)//TRIM(fext)
    ELSE IF(i<100) THEN
       WRITE(num2,fmt='(I2)') i
       number_file=TRIM(fbase)//'_'//TRIM(num2)//TRIM(fext)
    ELSE IF(i<1000) THEN
       WRITE(num3,fmt='(I3)') i
       number_file=TRIM(fbase)//'_'//TRIM(num3)//TRIM(fext)
    ELSE IF(i<10000) THEN
       WRITE(num4,fmt='(I4)') i
       number_file=TRIM(fbase)//'_'//TRIM(num4)//TRIM(fext)
    END IF

  END FUNCTION number_file

  SUBROUTINE read_input(m,x,v,n,file_name)

    IMPLICIT NONE
    REAL*16, ALLOCATABLE, INTENT(OUT) :: m(:), x(:,:), v(:,:)
    CHARACTER(len=256), INTENT(IN) :: file_name
    INTEGER, INTENT(OUT) :: n

    !Reads the input file

    n=file_length(file_name)

    WRITE(*,*) 'Particles in simulation:', n

    ALLOCATE(x(3,n),v(3,n),m(n))

    WRITE(*,*) 'Reading input'

    OPEN(7,file=file_name)
    DO i=1,n
       
       READ(7,*) m(i), x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i)

       !If only using 2 dimensions then ensure z compoents are set to zero
       IF(idim==2) THEN
          x(3,i)=0.d0
          v(3,i)=0.d0
       END IF
       
    END DO
    CLOSE(7)

    WRITE(*,*) 'Finished reading input'
    WRITE(*,*)

  END SUBROUTINE read_input

  FUNCTION cm(n,m,x)

    IMPLICIT NONE
    REAL*16, INTENT(IN) :: m(:), x(:,:)
    REAL*16 :: cm(3)
    INTEGER :: i, j, n

    !Calculate the CM vector from all particles
    
    cm=0.d0

    DO j=1,3
       DO i=1,n
          cm(j)=cm(j)+m(i)*x(j,i)
       END DO
    END DO

    cm=cm/SUM(m)

  END FUNCTION cm

  FUNCTION energy(m,x,v)

    IMPLICIT NONE
    INTEGER :: n, i, j
    REAL*16 :: kin, pot, energy, d
    REAL*16, INTENT(IN) :: m(:), x(:,:), v(:,:)

    !Calculates the total energy of all particles

    kin=0.d0
    pot=0.d0

    n=SIZE(x,2)

    !Calculate kinetic energy!
    DO i=1,n
       kin=kin+m(i)*((v(1,i)**2.d0)+(v(2,i)**2.d0)+(v(3,i)**2.d0))/2.d0
    END DO

!    WRITE(*,*) 'kin:', kin

    !Calculate potential energy
    DO i=1,n
       DO j=1,i-1
          IF(idim==2) THEN
             d=dist_2D(x(:,i),x(:,j))
          ELSE IF(idim==3) THEN
             d=dist_3D(x(:,i),x(:,j))
          END IF
          pot=pot+G*m(i)*m(j)*potential(d)
       END DO
    END DO

!    WRITE(*,*) 'pot', pot, d

    !Total is sum
    energy=kin+pot

  END FUNCTION energy

  FUNCTION potential(r)

    IMPLICIT NONE
    REAL*16 :: potential
    REAL*16, INTENT(IN) :: r
    REAL*16 :: rr

    !Potential energy without the G*m1*m2 pre-factor

    rr=r+soft

    IF(iforce==1) THEN
       potential=-1./rr
    ELSE IF(iforce==2) THEN
       potential=-1./rr**2
    ELSE IF(iforce==3) THEN
       potential=rr+rr**(-1)
    ELSE IF(iforce==4) THEN
       potential=rr**2.+rr**(-2)
    ELSE IF(iforce==5) THEN
       potential=rr**0.5+rr**(-0.5)
    ELSE IF(iforce==6) THEN
       potential=rr**(-3)-rr
    !ELSE IF(iforce==7) THEN
    !   potential=-1./(r+soft)
    ELSE
       STOP 'FORCE: Error, iforce not specified correctly'
    END IF
    
  END FUNCTION potential

  FUNCTION force(r)

    IMPLICIT NONE
    REAL*16 :: force
    REAL*16, INTENT(IN) :: r
    REAL*16 :: rr

    rr=r+soft

    !Force without the G*m1*m2 pre-factor
    !Force = -Grad V(r) (NB. *MINUS* grad V)
    !Assumes V is a function of r only

    IF(iforce==1) THEN
       force=-1./rr**2
    ELSE IF(iforce==2) THEN
       force=-2./rr**3
    ELSE IF(iforce==3) THEN
       force=-1.+rr**(-2)
    ELSE IF(iforce==4) THEN
       force=-2.*rr+2*rr**(-3)
    ELSE IF(iforce==5) THEN
       force=-0.5*rr**(-0.5)+0.5*rr**(-1.5)
    ELSE IF(iforce==6) THEN
       force=3.*rr**(-4)-rr**(-2)
    !ELSE IF(iforce==7) THEN
    !   potential=-1./(rr+soft)**2
    ELSE
       STOP 'FORCE: Error, iforce not specified correctly'
    END IF
    
  END FUNCTION force

  FUNCTION dist_2D(x1,x2)

    IMPLICIT NONE
    REAL*16 :: dist_2D, x1(2), x2(2)
    INTEGER :: i

    !Compute the distance between two 2D vectors

    dist_2D=0.d0

    DO i=1,2
       dist_2D=dist_2D+((x1(i)-x2(i))**2.d0)
    END DO

    dist_2D=sqrt(dist_2D)

  END FUNCTION dist_2D
  
  FUNCTION dist_3D(x1,x2)

    IMPLICIT NONE
    REAL*16 :: dist_3D, x1(3), x2(3)
    INTEGER :: i

    !Compute the distance between two 3D vectors

    dist_3D=0.d0

    DO i=1,3
       dist_3D=dist_3D+((x1(i)-x2(i))**2.d0)
    END DO

    dist_3D=sqrt(dist_3D)

  END FUNCTION dist_3D

  FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=256) :: file_name
    INTEGER ::n, file_length

    !Figures out the length of a file

    WRITE(*,*) 'File length of: ', TRIM(file_name)
    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*, end=301)
    END DO

301 CLOSE(7)

    n=n-1
    file_length=n

    WRITE(*,*) 'File length is:', file_length
    WRITE(*,*) ''

  END FUNCTION file_length

END PROGRAM nbody
