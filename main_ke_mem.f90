
! this program is used to data analysis for DNS data for .h5 files
! -Shaowu Pan. V2.0.0 - Jan.19
! -- new feature: it can read and write many files and one time!!!
! -- modify the the ke to limited in N/2 wave number

! -- command: h5fc -O3 fft.f90 fftshift.f90  main_ke.f95 -o fftke
! new new version h5fc -O3 fft.f90 fftshift.f90  main_ke_mem.f90 -o fftke

PROGRAM H5_MAIN_REFINEMENT
  use omp_lib
  USE HDF5 ! This module contains all necessary modules
  USE fft_mod

  IMPLICIT NONE

  CHARACTER(LEN=52) :: notmp,no

  CHARACTER(LEN=52), ALLOCATABLE :: filename_r(:)

  INTEGER(HID_T),ALLOCATABLE :: file_id_r(:)       ! File identifier
  INTEGER(HID_T),ALLOCATABLE :: dset_id_r(:,:)       ! Dataset identifier
  INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

  INTEGER     ::  N_REFINED,nmaxmax,nminmin
  INTEGER     ::  error ! Error flag
  INTEGER     ::  i, j, k,l,p1,k1,k2,k3
  INTEGER     ::  nmax,nmin
  INTEGER(kind=4),allocatable :: final_index(:)

  DOUBLE PRECISION, ALLOCATABLE :: data_in(:,:,:,:) ! Data buffers : because global.f90 has already defined it.
  DOUBLE PRECISION              :: dx,dy,dz,distance
! --- below: revise
  DOUBLE PRECISION, ALLOCATABLE :: final(:),output(:)
  DOUBLE PRECISION              :: su_output,macht,urms2,k0

  COMPLEX(kind=dp), ALLOCATABLE :: fu(:,:,:),fv(:,:,:),fw(:,:,:)
  REAL(kind=dp), ALLOCATABLE :: afu(:,:,:),afv(:,:,:),afw(:,:,:)
 ! COMPLEX, ALLOCATABLE :: fk1(:,:,:,:),fk2(:,:,:,:),fk3(:,:,:,:)
!-------------------- I: INPUT: RAW DATA INFORMATION -------------------------------------------
  INTEGER           ::  N_RAW    = 11              ! ------------no need--modify
  INTEGER           ::  N_MESH   = 64             ! grid size
  CHARACTER(LEN=3) :: dsetname(6) ! Dataset Name

!-------------------- II: INPUT: Refinement degree ---------------------------------------------
  INTEGER           ::  L_REFINE = 0           ! --------------modify

!-------------------- III: ALLOCALATION AND INITILAIZATION---------------------------------------
  do i=1,3;
      data_dims(i) = N_MESH       ! Data_Dimension
  end do
  dx=2d0*pi/data_dims(1); dy=2d0*pi/data_dims(2); dz=2d0*pi/data_dims(3);

  DSETNAME(1)='den' ! density
  DSETNAME(2)='xve' ! U should be zero
  DSETNAME(3)='pre' !
  DSETNAME(4)='mas'
  DSETNAME(5)='yve'
  DSETNAME(6)='zve'

  call system('mkdir ke_spectrum')

  N_REFINED = L_REFINE*(N_RAW - 1) + 1  ! formula for 21
!
  ALLOCATE(filename_r(N_RAW))
  ALLOCATE(file_id_r(N_RAW))
  ALLOCATE(dset_id_r(N_RAW,6))


  k0=4d0
  macht=0.1d0
  urms2=11.2d0*pi**2*macht**2
  k0=4;
!#

!  -------------------- IV: NAME of FOLDER---------------------------------------


  do l=1,N_RAW
      write(notmp,*) l-1
    IF((l-1) < 10)THEN
      no = '00'//TRIM(ADJUSTL(notmp))
    ELSEIF((l-1) >= 100)THEN
      no = TRIM(ADJUSTL(notmp))
    ELSE
      no = '0'//TRIM(ADJUSTL(notmp))
    END IF
    filename_r(l) = "./output/outputhdf"//TRIM(ADJUSTL(no))//".h5"
  end do


!--- Create do loops for each file
!----------------------data analysis--beigns-  --------------------------------

! -- init for HDF lib


do l=1,N_RAW
! --- allocate data_in. = 4GB for 512^3 is used here
call h5open_f(error);
   ALLOCATE(data_in(data_dims(1),data_dims(2),data_dims(3),6))

! -- 9G

!---------------------- box 1: READING DATA BEGIN---------------------------------------------------
   CALL h5fopen_f (filename_r(l), H5F_ACC_RDWR_F, file_id_r(l), error)   ! 1. in: filename out: file_id

      do i=1,6; ! loop for each component
              !--I. Open an existing dataset. ! must tell the code which dataset.
          CALL h5dopen_f(file_id_r(l), dsetname(i), dset_id_r(l,i), error) ! 2. in: file_id,dsetname out:dset_id
          !
          CALL h5dread_f(dset_id_r(l,i), H5T_NATIVE_DOUBLE, data_in(:,:,:,i), data_dims, error)! Read the dataset.

          CALL h5dclose_f(dset_id_r(l,i), error) !end -1

      end do ! loop for each component -  end

   CALL h5fclose_f(file_id_r(l), error)
CALL h5close_f(error) ! close the whole HDF5 lib
!--------------------- box 1: READING DATA END----------------------------------------------------

!   CALL ens_cal(kd,ens1,ens2,ens3,data_dims) ! the ens is calcualted with 0.5
!deallocate(Data_in); 
! -- 3G
ALLOCATE(fu(data_dims(1),data_dims(2),data_dims(3)))  

!$OMP PARALLEL WORKSHARE
   fu(:,:,:)=DCMPLX(dsqrt(data_in(:,:,:,1))*data_in(:,:,:,2),0);
!$OMP END PARALLEL WORKSHARE

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fu(:,j,k));
       fu(:,j,k)=fu(:,j,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fu(j,:,k));
       fu(j,:,k)=fu(j,:,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fu(j,k,:));
       fu(j,k,:)=fu(j,k,:)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

allocate(afu(data_dims(1),data_dims(2),data_dims(3)));
!$OMP PARALLEL WORKSHARE
afu=abs(fu);
!$OMP END PARALLEL WORKSHARE
deallocate(fu);

!=-------============ 2- ens2
   
ALLOCATE(fv(data_dims(1),data_dims(2),data_dims(3)))
!$OMP PARALLEL WORKSHARE
    fv(:,:,:)=DCMPLX(dsqrt(data_in(:,:,:,1))*data_in(:,:,:,5),0);
!$OMP END PARALLEL WORKSHARE

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fv(:,j,k));
       fv(:,j,k)=fv(:,j,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fv(j,:,k));
       fv(j,:,k)=fv(j,:,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fv(j,k,:));
       fv(j,k,:)=fv(j,k,:)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

allocate(afv(data_dims(1),data_dims(2),data_dims(3)));
!$OMP PARALLEL WORKSHARE
afv=abs(fv);
!$OMP END PARALLEL WORKSHARE
deallocate(fv);
    ! -- calcualte k's module
! =================== 3- ens3
    
ALLOCATE(fw(data_dims(1),data_dims(2),data_dims(3))) 
!$OMP PARALLEL WORKSHARE
    fw(:,:,:)=DCMPLX(dsqrt(data_in(:,:,:,1))*data_in(:,:,:,6),0);
!$OMP END PARALLEL WORKSHARE
deallocate(Data_in)

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fw(:,j,k));
       fw(:,j,k)=fw(:,j,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fw(j,:,k));
       fw(j,:,k)=fw(j,:,k)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(fw(j,k,:));
       fw(j,k,:)=fw(j,k,:)/N_MESH
    end do;end do
!$OMP END PARALLEL DO

allocate(afw(data_dims(1),data_dims(2),data_dims(3)));
!$OMP PARALLEL WORKSHARE
afw=abs(fw);
!$OMP END PARALLEL WORKSHARE
deallocate(fw);

!---------- calcualted all the necessary things

  !-- conversion should happen before
  
  ! -- at there for 512^3 case, it should occupy memory
  !
  !  3G
  !  
  ! -- afu,afv,afw

 
  ALLOCATE(final(N_MESH**3))
  allocate(final_index(N_mesh**3))
  final=0;
  final_index=0; 
  ! 
  ! 5G
  !

    p1=0;

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,k1,k2,k3,distance) LASTPRIVATE(p1)
    do i=1,N_MESH;do j=1,N_MESH;do k=1,N_MESH;

       p1 = (i-1)*N_MESH**2 + (j-1)*N_MESH + k

       ! - box-1. special treatment for doing specturm operations for FFT.
       if (i.lt.(N_MESH/2+1))then
           k1=i-1
       elseif (i.eq.N_MESH/2+1) then
           k1=0
       else
           k1=i-1-N_MESH
       endif

       if (j.lt.(N_MESH/2+1))then
           k2=j-1
       elseif (j.eq.N_MESH/2+1)then
           k2=0
       else
           k2=j-1-N_MESH
       endif

       if (k.lt.(N_MESH/2+1))then
           k3=k-1
       elseif (k.eq.N_MESH/2+1)then
           k3=0
       else
           k3=k-1-N_MESH
       endif
       ! =============distance =====================
       distance=(k1*k1)+(k2*k2)+(k3*k3);
       ! -----============put them into =======================
       final_index(p1) = idnint(dsqrt(distance));
                         ! nint(dsqrt(dble(i-(N_MESH/2+1))**2+ &
                         ! dble(j-(N_MESH/2+1))**2+dble(k-(N_MESH/2+1))**2))
       if (((k1.eq.0).AND.(k2.eq.0)).AND.(k3.eq.0)) THEN
          final(p1)=0d0
       else
          final(p1)=(afu(i,j,k)**2+afv(i,j,k)**2+afw(i,j,k)**2)/2d0 !4d0*pi*distance/2d0
          !final(p1,l)=final(p1,l)/N_MESH**3
       endif
    enddo;enddo;enddo
!$OMP END DO
!$OMP END PARALLEL

deallocate(afu)
deallocate(afv)
deallocate(afw)
! --2G



! ----- part two: data output 

!$OMP PARALLEL WORKSHARE
nmax=maxval(final_index(:));
nmin=minval(final_index(:));
!$OMP END PARALLEL WORKSHARE

nmaxmax=nmax;
nminmin=nmin;

allocate(output(nminmin:nmaxmax));

output=0d0;

! -- 2G at most

!$OMP PARALLEL DO PRIVATE(p1) 
  do p1=1,N_MESH**3
!$OMP ATOMIC     
     output(final_index(p1)) = output(final_index(p1)) + final(p1);
  end do
!$OMP END PARALLEL DO

deallocate(final);
deallocate(final_index);


  su_output=0d0
  do i=nminmin,nmaxmax
     su_output=su_output+output(i);
  enddo;
  
  write(notmp,*) l-1
  IF((l-1) < 10)THEN
      no = '00'//TRIM(ADJUSTL(notmp))
  ELSEIF((l-1) >= 100)THEN
      no = TRIM(ADJUSTL(notmp))
  ELSE
      no = '0'//TRIM(ADJUSTL(notmp))
  END IF

! 123 FORMAT(2E18.7)
  open(unit=20,file='./ke_spectrum/Lab'//TRIM(ADJUSTL(no))//".dat")!
  write(20,*) "VARIABLES='wave number','Ek'"

  do i=nminmin,nmaxmax
     if (i.lt.(N_mesh*sqrt(1d0)/2))then
     if (output(i).gt.0.000000000000d0) then
         write(20,*) i,output(i),urms2/k0**5*16d0*dsqrt(2d0/pi)*i**4*exp(-2d0*i**2/k0**2)
     endif;end if
  end do
  close(20)

deallocate(output)
!
!-0G
!

end do




END PROGRAM H5_MAIN_REFINEMENT


