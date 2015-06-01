PROGRAM nbody

  IMPLICIT NONE
  REAL*16, ALLOCATABLE :: xi(:,:), vi(:,:), m(:), x(:,:), v(:,:), xnew(:,:), vnew(:,:), xres(:,:,:), tres(:), vres(:,:,:)
  INTEGER :: i, j, k, np, icm
  INTEGER :: iE, IL
  REAL*16 :: tf, dt, for(3), e_init, e_final, xcm(3), vcm(3), l_init(3), l_final(3), e, enew, de, t, dtmax, eratio
  REAL*16 :: L(3), Lnew(3), Lmod, Lmodnew, Lratio, dL, acc
  REAL*16, PARAMETER :: G=1.d0, pi=3.141592654d0, soft=0.d0, ti=0.d0!, acc=1.e-6
  INTEGER, PARAMETER :: n=1001
  CHARACTER(len=64) :: infile, time, boost, accuracy
  LOGICAL :: mr_logic

  CALL get_command_argument(1,infile)
  IF(infile=='') STOP 'You need to specify an input file'

  CALL get_command_argument(2,time)
  IF(time=='') STOP 'You need to specify a simulation time length (yrs)'
  READ(time,*) tf

  CALL get_command_argument(3,boost)
  IF(boost=='') STOP 'You need to specify CM boost (1) or not (0)'
  READ(boost,*) icm
  IF((icm .NE. 1) .AND. (icm .NE. 0)) STOP 'You need to specify CM boost (1) or not (0)'

  CALL get_command_argument(4,accuracy)
  IF(accuracy=='') STOP 'You need to specify an accuracy parameter (e.g. 1e-6)'
  READ(accuracy,*) acc

  INQUIRE(file=infile,exist=mr_logic)
  IF(mr_logic .EQV. .FALSE.) STOP 'Input file does not exist'

  WRITE(*,*)
  WRITE(*,*) 'N-body integration code'
  WRITE(*,*) '======================='
  WRITE(*,*)

  CALL read_input(m,xi,vi,np,infile)

  xcm=cm(np,m,xi)
  vcm=cm(np,m,vi)

!  WRITE(*,*) 'Move to CM position?'
!  WRITE(*,*) '1 - Yes'
!  WRITE(*,*) '0 - No'
!  READ(*,*) icm
!  WRITE(*,*)
  WRITE(*,*) 'Boosting to CM position'

  IF(icm==1) THEN
     DO i=1,np
        DO j=1,3
           xi(j,i)=xi(j,i)-xcm(j)
           vi(j,i)=vi(j,i)-vcm(j)
        END DO
     END DO
  END IF

  WRITE(*,*) 'CM move complete'
  WRITE(*,*)

  WRITE(*,*) 'Mass in solar masses'
  WRITE(*,*) 'Length units in AU'
  WRITE(*,*) 'Velocity in 2*pi*AU/year'
!  WRITE(*,*)
  WRITE(*,*) 'Time of simulation [years]', tf
!  READ(*,*) tf
  tf=tf*2.d0*pi
  WRITE(*,*)

  ALLOCATE(x(3,np),xnew(3,np),v(3,np),vnew(3,np),xres(3,np,n),vres(3,np,n),tres(n))

  x=0.d0
  xnew=0.d0
  v=0.d0
  vnew=0.d0
  t=ti

  DO k=1,3
     DO i=1,np
        x(k,i)=xi(k,i)
        v(k,i)=vi(k,i)
        xres(k,i,1)=x(k,i)
        vres(k,i,1)=v(k,i)
     END DO
  END DO

  iL=1
  iE=1

  L=am(m,x,v)
  Lmod=length(L)
  IF(Lmod<=.0000001d0) iL=0

  dt=(tf-ti)/float(n-1)
  dtmax=dt

  WRITE(*,*) 'Starting intergation'
  IF(iE==1) WRITE(*,*) 'Integrating while checking energy conservation'
  IF(iL==1) WRITE(*,*) 'Integrating while checking angular momentum conservation'
  WRITE(*,*) 'Accuracy parameter:', acc
  WRITE(*,*)

  i=1

  DO

     IF(iE==1) e=energy(m,x,v)
     IF(iL==1) THEN
        L=am(m,x,v)
        Lmod=length(L)
     END IF

     CALL rk4(np,x(:,:),xnew(:,:),v(:,:),vnew(:,:),dt)

     IF(iE==1) enew=energy(m,xnew,vnew)
     IF(iL==1) THEN
        Lnew=am(m,x,v)
        Lmodnew=length(Lnew)
     END IF

     IF(iE==1) eratio=ABS(-1.d0+enew/e)
     IF(iE==1) de=acc*dt/tf
     IF(iL==1) Lratio=ABS(-1.d0+Lmodnew/Lmod)
     IF(iL==1) dL=acc*dt/tf

     !WRITE(7,fmt='(4F20.10)') t, dt, eratio, de
     !WRITE(*,fmt='(4F20.10)') t, dt, eratio, de
     IF((iE==0 .OR. eratio<de) .AND. (iL==0 .OR. Lratio<dL)) THEN
        x=xnew
        v=vnew
        t=t+dt
        IF(t>=dtmax*float(i)) THEN
           IF(float(i)>float(n)/10. .AND. float(i)<float(n)/10.+1.) WRITE(*,*) '10% Complete'
           IF(float(i)>2.*float(n)/10. .AND. float(i)<2.*float(n)/10.+1.) WRITE(*,*) '20% Complete'
           IF(float(i)>3.*float(n)/10. .AND. float(i)<3.*float(n)/10.+1.) WRITE(*,*) '30% Complete'
           IF(float(i)>4.*float(n)/10. .AND. float(i)<4.*float(n)/10.+1.) WRITE(*,*) '40% Complete'
           IF(float(i)>5.*float(n)/10. .AND. float(i)<5.*float(n)/10.+1.) WRITE(*,*) '50% Complete'
           IF(float(i)>6.*float(n)/10. .AND. float(i)<6.*float(n)/10.+1.) WRITE(*,*) '60% Complete'
           IF(float(i)>7.*float(n)/10. .AND. float(i)<7.*float(n)/10.+1.) WRITE(*,*) '70% Complete'
           IF(float(i)>8.*float(n)/10. .AND. float(i)<8.*float(n)/10.+1.) WRITE(*,*) '80% Complete'
           IF(float(i)>9.*float(n)/10. .AND. float(i)<9.*float(n)/10.+1.) WRITE(*,*) '90% Complete'
           i=i+1
           xres(:,:,i)=x(:,:)
           vres(:,:,i)=v(:,:)
           tres(i)=t
           IF(i==n) THEN
              WRITE(*,*) '100% Complete'
              EXIT
           END IF
        END IF
        IF(dt<dtmax/1.1d0) THEN
           dt=dt*2.d0
        END IF
     ELSE

        dt=dt/2.d0

        IF(dt<dtmax/1.e16) THEN
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

  WRITE(*,*) 'Initial energy:', e_init
  WRITE(*,*) 'Final energy:', e_final
  WRITE(*,*) 'Ratio:', e_final/e_init
  WRITE(*,*)

  CALL results(np,1000,tres,xres,vres,m)

CONTAINS

  FUNCTION am(m,x,v)

    IMPLICIT NONE
    REAL*16 :: am(3)
    REAL*16, INTENT(IN) :: m(:), x(:,:), v(:,:)
    INTEGER :: i, n

    am=0.d0

    n=SIZE(x,2)

    DO i=1,n
       !       am(1)=am(1)+m(i)*(x(2,i)*v(3,i)-x(3,i)*v(2,i))
       !       am(2)=am(2)+m(i)*(x(3,i)*v(1,i)-x(1,i)*v(3,i))
       !      am(3)=am(3)+m(i)*(x(1,i)*v(2,i)-x(2,i)*v(1,i))
       am=am+m(i)*cross_product(x(:,i),v(:,i))
    END DO

  END FUNCTION am

  SUBROUTINE matrix(n,x,dij)

    IMPLICIT NONE
    REAL*16, INTENT(IN) :: x(:,:)
    INTEGER :: n
    REAL*16, INTENT(OUT) :: dij(:,:,:)
    INTEGER :: i, j, k
    REAL*16 :: d3

    dij=0.d0

    DO i=1,n
       DO j=1,i-1
          d3=(dist(x(:,i),x(:,j))+soft)**3.d0
          DO k=1,3
             dij(k,i,j)=-(x(k,i)-x(k,j))/d3
             dij(k,j,i)=-dij(k,i,j)
          END DO
       END DO
    END DO

    !Hopefully the below part isn't necessary!

    !    DO i=1,n
    !       DO k=1,3
    !          dij(k,i,i)=0.d0
    !       END DO
    !    END DO

    !    dij=dij*G

  END SUBROUTINE matrix

  SUBROUTINE rk4(n,xin,xout,vin,vout,dt)

    REAL*16, INTENT(IN) :: xin(:,:), vin(:,:), dt
    REAL*16, INTENT(OUT) :: xout(:,:), vout(:,:)
    INTEGER, INTENT(IN) :: n
    REAL*16, ALLOCATABLE :: kx1(:,:), kx2(:,:), kx3(:,:), kx4(:,:), kv1(:,:), kv2(:,:), kv3(:,:), kv4(:,:), disp(:,:,:)
    REAL*16 :: accl
    INTEGER :: i, j, k

    ALLOCATE(kx1(3,n),kx2(3,n),kx3(3,n),kx4(3,n),kv1(3,n),kv2(3,n),kv3(3,n),kv4(3,n))
    ALLOCATE(disp(3,n,n))

    kv1=0.d0
    kv2=0.d0
    kv3=0.d0
    kv4=0.d0

    kx1=0.d0
    kx2=0.d0
    kx3=0.d0
    kx4=0.d0

    CALL matrix(n,xin,disp)

    DO k=1,3
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=m(j)*disp(k,i,j)
                kv1(k,i)=kv1(k,i)+accl*dt
             END IF
          END DO
          kx1(k,i)=vin(k,i)*dt
       END DO
    END DO

    CALL matrix(n,xin+kx1/2.d0,disp)
    DO k=1,3
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=m(j)*disp(k,i,j)
                kv2(k,i)=kv2(k,i)+accl*dt
             END IF
          END DO
          kx2(k,i)=(vin(k,i)+kv1(k,i)/2.d0)*dt
       END DO
    END DO

    CALL matrix(n,xin+kx2/2.d0,disp)
    DO k=1,3
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=m(j)*disp(k,i,j)
                kv3(k,i)=kv3(k,i)+accl*dt
             END IF
          END DO
          kx3(k,i)=(vin(k,i)+kv2(k,i)/2.d0)*dt
       END DO
    END DO

    CALL matrix(n,xin+kx3,disp)
    DO k=1,3
       DO i=1,n
          DO j=1,n
             IF(i .NE. j) THEN
                accl=m(j)*disp(k,i,j)
                kv4(k,i)=kv4(k,i)+accl*dt
             END IF
          END DO
          kx4(k,i)=(vin(k,i)+kv3(k,i))*dt
       END DO
    END DO

    DEALLOCATE(disp)

    DO k=1,3
       DO i=1,n
          xout(k,i)=xin(k,i)+(kx1(k,i)+(2.d0*kx2(k,i))+(2.d0*kx3(k,i))+kx4(k,i))/6.d0
          vout(k,i)=vin(k,i)+(kv1(k,i)+(2.d0*kv2(k,i))+(2.d0*kv3(k,i))+kv4(k,i))/6.d0
       END DO
    END DO

    DEALLOCATE(kx1,kx2,kx3,kx4,kv1,kv2,kv3,kv4)

  END SUBROUTINE rk4

  FUNCTION length(x)

    IMPLICIT NONE
    REAL*16 :: length
    REAL*16, INTENT(IN) :: x(3)

    length=sqrt(x(1)**2.+x(2)**2.+x(3)**2.)

  END FUNCTION length

  FUNCTION dot_product(x,y)

    IMPLICIT NONE
    REAL*16 :: dot_product
    REAL*16, INTENT(IN) :: x(3), y(3)

    dot_product=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)

  END FUNCTION dot_product

  FUNCTION cross_product(x,y)

    IMPLICIT NONE
    REAL*16 :: cross_product(3)
    REAL*16, INTENT(IN) :: x(3), y(3)

    cross_product(1)=x(2)*y(3)-x(3)*y(2)
    cross_product(2)=x(3)*y(1)-x(1)*y(3)
    cross_product(3)=x(1)*y(2)-x(2)*y(1)

  END FUNCTION cross_product

  SUBROUTINE results(np,n,t,x,v,m)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, np
    REAL*16, ALLOCATABLE :: t(:), x(:,:,:), v(:,:,:), m(:)
    CHARACTER(len=1) :: file_num1, part_num1
    CHARACTER(len=2) :: file_num2, part_num2
    CHARACTER(len=3) :: file_num3, part_num3
    CHARACTER(len=4) :: file_num4, part_num4
    CHARACTER(len=64) :: fname, stem, ext
    INTEGER :: tail, i, j, k

    t=t/(2.d0*pi)

    stem='./data/particle'
    ext='.dat'

    DO j=1,np
!       WRITE(file_num1,fmt='(I1)') j
       fname=number_file(stem,j,ext)
 !      fname='./data/particle_'//part_num//'.dat'
       OPEN(7,file=fname)
       DO i=1,n
          WRITE(7,fmt='(7F20.10)') t(i), x(1,j,i), x(2,j,i), x(3,j,i), v(1,j,i), v(2,j,i), v(3,j,i)
       END DO
       CLOSE(7)
    END DO

    OPEN(7,file='end.dat')
    DO i=1,np
       WRITE(7,fmt='(7F15.7)') m(i), x(1,i,n), x(2,i,n), x(3,i,n), v(1,i,n), v(2,i,n), v(3,i,n)
    END DO
    CLOSE(7)

    !Tail governs the tail length of particles in the dump
!    tail=1

!    stem='./data/dump'
!    ext=''

!    DO k=1,np
!       DO i=1,n

!          WRITE(part_num,fmt='(I1)') k

 !         IF(i<10) THEN
 !            WRITE(file_num1,fmt='(I1)') i
 !            fname='./data/dump_'//part_num//'_'//file_num1
 !         ELSE IF(i<100) THEN
 !            WRITE(file_num2,fmt='(I2)') i
 !            fname='./data/dump_'//part_num//'_'//file_num2
 !         ELSE IF(i<1000) THEN
 !            WRITE(file_num3,fmt='(I3)') i
 !            fname='./data/dump_'//part_num//'_'//file_num3
 !         ELSE IF(i<10000) THEN
 !            WRITE(file_num4,fmt='(I4)') i
 !            fname='./data/dump_'//part_num//'_'//file_num4
 !         END IF

!          fname=number_file2(stem,k,i,ext)

!          OPEN(7,file=fname)
!          DO j=i,i-tail+1,-1
!             IF(j<1) EXIT
!             WRITE(7,fmt='(4F20.10)') t(j), x(1,k,j), x(2,k,j), x(3,k,j)
!          END DO
!          CLOSE(7)

!       END DO
    ! END DO

  END SUBROUTINE results

   FUNCTION number_file(fbase,i,fext)

    IMPLICIT NONE
    CHARACTER(len=64) number_file, fbase, fext
    CHARACTER(len=4) num4
    CHARACTER(len=3) num3
    CHARACTER(len=2) num2
    CHARACTER(len=1) num1
    INTEGER :: i

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

  FUNCTION number_file2(fbase,i,j,fext)

    IMPLICIT NONE
    CHARACTER(len=64) number_file2, fbase, fext
    CHARACTER(len=4) num4, rum4
    CHARACTER(len=3) num3, rum3
    CHARACTER(len=2) num2, rum2
    CHARACTER(len=1) num1, rum1
    INTEGER :: i, j

    IF(i<10) THEN

       WRITE(num1,fmt='(I1)') i

       IF(j<10) THEN
          WRITE(rum1,fmt='(I1)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num1)//'_'//TRIM(rum1)//TRIM(fext)
       ELSE IF(j<100) THEN
          WRITE(rum2,fmt='(I2)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num1)//'_'//TRIM(rum2)//TRIM(fext)
       ELSE IF(j<1000) THEN
          WRITE(rum3,fmt='(I3)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num1)//'_'//TRIM(rum3)//TRIM(fext)
       ELSE IF(j<10000) THEN
          WRITE(rum4,fmt='(I4)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num1)//'_'//TRIM(rum4)//TRIM(fext)
       END IF

    ELSE IF(i<100) THEN

       WRITE(num2,fmt='(I2)') i

       IF(j<10) THEN
          WRITE(rum1,fmt='(I1)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num2)//'_'//TRIM(rum1)//TRIM(fext)
       ELSE IF(j<100) THEN
          WRITE(rum2,fmt='(I2)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num2)//'_'//TRIM(rum2)//TRIM(fext)
       ELSE IF(j<1000) THEN
          WRITE(rum3,fmt='(I3)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num2)//'_'//TRIM(rum3)//TRIM(fext)
       ELSE IF(j<10000) THEN
          WRITE(rum4,fmt='(I4)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num2)//'_'//TRIM(rum4)//TRIM(fext)
       END IF
  
    ELSE IF(i<1000) THEN

       WRITE(num3,fmt='(I3)') i

       IF(j<10) THEN
          WRITE(rum1,fmt='(I1)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num3)//'_'//TRIM(rum1)//TRIM(fext)
       ELSE IF(j<100) THEN
          WRITE(rum2,fmt='(I2)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num3)//'_'//TRIM(rum2)//TRIM(fext)
       ELSE IF(j<1000) THEN
          WRITE(rum3,fmt='(I3)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num3)//'_'//TRIM(rum3)//TRIM(fext)
       ELSE IF(j<10000) THEN
          WRITE(rum4,fmt='(I4)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num3)//'_'//TRIM(rum4)//TRIM(fext)
       END IF

    ELSE IF(i<10000) THEN

       WRITE(num4,fmt='(I4)') i

       IF(j<10) THEN
          WRITE(rum1,fmt='(I1)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num4)//'_'//TRIM(rum1)//TRIM(fext)
       ELSE IF(j<100) THEN
          WRITE(rum2,fmt='(I2)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num4)//'_'//TRIM(rum2)//TRIM(fext)
       ELSE IF(j<1000) THEN
          WRITE(rum3,fmt='(I3)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num4)//'_'//TRIM(rum3)//TRIM(fext)
       ELSE IF(j<10000) THEN
          WRITE(rum4,fmt='(I4)') j
          number_file2=TRIM(fbase)//'_'//TRIM(num4)//'_'//TRIM(rum4)//TRIM(fext)
       END IF

    END IF

  END FUNCTION number_file2

  SUBROUTINE read_input(m,x,v,n,file_name)

    REAL*16, ALLOCATABLE, INTENT(OUT) :: m(:), x(:,:), v(:,:)
    CHARACTER(len=64), INTENT(IN) :: file_name
    INTEGER, INTENT(OUT) :: n
    INTEGER :: inc

    n=file_length(file_name)

!    IF(file_name=='end.dat') n=n-1

    WRITE(*,*) 'Particles in simulation:', n

    ALLOCATE(x(3,n),v(3,n),m(n))

    WRITE(*,*) 'Reading input'

    OPEN(7,file=file_name)
    DO i=1,n
       READ(7,*) m(i), x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i)
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
    INTEGER :: n, i, j, k
    REAL*16 :: kin, pot, energy
    REAL*16 :: m(:), x(:,:), v(:,:)

    kin=0.d0
    pot=0.d0

    n=SIZE(x,2)

    !Calculate kinetic energy!

    DO i=1,n
       kin=kin+m(i)*((v(1,i)**2.d0)+(v(2,i)**2.d0)+(v(3,i)**2.d0))/2.d0
    END DO

    !Calculate potential energy

    DO i=1,n
       DO j=1,i-1
          pot=pot-G*m(i)*m(j)/(dist(x(:,i),x(:,j))+soft)
       END DO
    END DO

    energy=kin+pot

  END FUNCTION energy

  FUNCTION dist(x1,x2)

    IMPLICIT NONE
    REAL*16 :: dist, x1(3), x2(3)
    INTEGER :: i

    dist=0.d0

    DO i=1,3
       dist=dist+((x1(i)-x2(i))**2.d0)
    END DO

    dist=sqrt(dist)

  END FUNCTION dist

  FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=64) :: file_name
    INTEGER ::n, file_length
    REAL :: data

    WRITE(*,*) 'File length of:', file_name
    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*, end=301) data
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)

    n=n-1
    file_length=n

    WRITE(*,*) 'File length is:', file_length
    WRITE(*,*) ''

  END FUNCTION file_length

END PROGRAM nbody
