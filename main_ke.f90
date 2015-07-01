


! this program is used to data analysis for DNS data for .h5 files
! -Shaowu Pan. V2.0.0 - Jan.19
! -- new feature: it can read and write many files and one time!!!
! -- modify the the ke to limited in N/2 wave number

! -- command: h5fc -O3 fft.f90 fftshift.f90  main_ke.f95 -o fftke
! new new version h5fc -O3 fft.f90 fftshift.f90  main_ke.f90 -o fftke

PROGRAM H5_MAIN_FFT

  USE HDF5 ! This module contains all necessary modules
  USE fft_mod
  
  IMPLICIT NONE

  CHARACTER(LEN=52) :: notmp,no

  CHARACTER(LEN=52), ALLOCATABLE :: filename_r(:),filename_w(:)

  INTEGER(HID_T),ALLOCATABLE :: file_id_r(:)       ! File identifier
  INTEGER(HID_T),ALLOCATABLE :: file_id_w(:)       ! File identifier
  INTEGER(HID_T),ALLOCATABLE :: dset_id_r(:,:)       ! Dataset identifier
  INTEGER(HID_T),ALLOCATABLE :: dset_id_w(:,:)       ! Dataset identifier
  INTEGER(HID_T),ALLOCATABLE :: dspace_id_w(:)
  INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

  INTEGER     ::  N_REFINED,nmaxmax,nminmin
  INTEGER     ::  error ! Error flag
  INTEGER     ::  i, j, k,l,p1,p2,k1,k2,k3,rank=3
  INTEGER(kind=8),allocatable :: nmax(:),nmin(:),final_index(:,:),output_index(:,:)
  integer     :: final_index_1,final_index_2,psum  

  DOUBLE PRECISION, ALLOCATABLE :: buff(:,:,:,:) !, ALLOCATABLE :: buff(:,:,:) ! Data buffers
  DOUBLE PRECISION, ALLOCATABLE :: data_in(:,:,:,:,:) ! Data buffers
  DOUBLE PRECISION, ALLOCATABLE :: data_out(:,:,:,:,:) ! Data buffers
  DOUBLE PRECISION, ALLOCATABLE :: final(:,:),output(:,:),cum_output(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: output_value(:,:),su_output(:)

  double precision :: distance,macht,urms2,k0
  complex(kind=dp), ALLOCATABLE :: f(:,:,:),fu(:,:,:),fv(:,:,:),fw(:,:,:)
 ! complex(kind=dp), ALLOCATABLE :: fk1(:,:,:,:),fk2(:,:,:,:),fk3(:,:,:,:)
!-------------------- I: INPUT: RAW DATA INFORMATION -------------------------------------------
  INTEGER           ::  N_RAW    = 11              ! ------------no need--modify
  INTEGER           ::  N_MESH   = 64             ! grid size
  CHARACTER(LEN=3)  ::  dsetname(6) ! Dataset Name

!-------------------- II: INPUT: Refinement degree ---------------------------------------------
  INTEGER           ::  L_REFINE = 0           ! --------------modify

!-------------------- III: ALLOCALATION AND INITILAIZATION---------------------------------------
  do i=1,3;       
      data_dims(i) = N_MESH       ! Data_Dimension
  end do

  call system('mkdir ke_spectrum')

  DSETNAME(1)='den' ! density
  DSETNAME(2)='xve' ! U should be zero
  DSETNAME(3)='pre' !
  DSETNAME(4)='mas'
  DSETNAME(5)='yve'
  DSETNAME(6)='zve'

  N_REFINED = L_REFINE*(N_RAW - 1) + 1  ! formula for 21
! 
  ALLOCATE(filename_r(N_RAW))
  ALLOCATE(filename_w(N_REFINED))
  ALLOCATE(file_id_r(N_RAW))
  ALLOCATE(file_id_w(N_REFINED))
  ALLOCATE(dset_id_r(N_RAW,6))
  ALLOCATE(dset_id_w(N_REFINED,6))
  ALLOCATE(dspace_id_w(N_REFINED))
  ALLOCATE(buff(data_dims(1),data_dims(2),data_dims(3),6))
  ALLOCATE(data_in(N_RAW,data_dims(1),data_dims(2),data_dims(3),6))
  ALLOCATE(data_out(N_REFINED,data_dims(1),data_dims(2),data_dims(3),6))
  ALLOCATE(f(data_dims(1),data_dims(2),data_dims(3)))
  ALLOCATE(fu(data_dims(1),data_dims(2),data_dims(3)))
  ALLOCATE(fv(data_dims(1),data_dims(2),data_dims(3)))
  ALLOCATE(fw(data_dims(1),data_dims(2),data_dims(3)))
  allocate(nmax(N_raw),nmin(n_raw))
  ALLOCATE(final(N_MESH**3,N_RAW))
  allocate(final_index(N_mesh**3,n_raw))
  allocate(su_output(1:N_RAW))
!-------------------- IV: NAME of FOLDER---------------------------------------


  do l=1,N_RAW
      write(notmp,*) l-1
    IF((l-1) < 10)THEN
      no = '00'//TRIM(ADJUSTL(notmp))
    ELSEIF((l-1) >= 100)THEN
      no = TRIM(ADJUSTL(notmp))
    ELSE
      no = '0'//TRIM(ADJUSTL(notmp))
    END IF
    filename_r(l) = "outputhdf"//TRIM(ADJUSTL(no))//".h5"
    !write(*,*)filename_r(l)
    su_output(l)=0;
  end do

  ! --------exact spectrum parameters
  macht=0.1d0
  urms2=11.2d0*pi**2*macht**2
  k0=4;
  ! ---------------------------------------------------------------------------------
  do l=1,N_REFINED

    write(notmp,*) l-1

    IF((l-1) < 10)THEN
      no = '00'//TRIM(ADJUSTL(notmp))
    ELSEIF((l-1) >= 100)THEN
      no = TRIM(ADJUSTL(notmp))
    ELSE
      no = '0'//TRIM(ADJUSTL(notmp))
    END IF
    filename_w(l) = "./refined/outputhdf"//TRIM(ADJUSTL(no))//".h5"

  end do

!---------------------- V: DATA ANALYSIS ---------------------------------------
  CALL h5open_f(error) !0. initialize the fortran interface

  do l=1,N_RAW;
!---------------------- box 1: READING DATA BEGIN---------------------------------------------------
  CALL h5fopen_f (filename_r(l), H5F_ACC_RDWR_F, file_id_r(l), error)   ! 1. in: filename out: file_id
  
      do i=1,6; ! loop for each component
              !--I. Open an existing dataset. ! must tell the code which dataset.
          CALL h5dopen_f(file_id_r(l), dsetname(i), dset_id_r(l,i), error) ! 2. in: file_id,dsetname out:dset_id
          !
          CALL h5dread_f(dset_id_r(l,i), H5T_NATIVE_DOUBLE, buff(:,:,:,i), data_dims, error)! Read the dataset.

	  data_in(l,:,:,:,i)=buff(:,:,:,i)

          CALL h5dclose_f(dset_id_r(l,i), error) !end -1

          !--I. end -----------
      end do ! loop for each component -  end

    CALL h5fclose_f(file_id_r(l), error)

!--------------------- box 1: READING DATA END----------------------------------------------------
  end do;
  
  CALL h5close_f(error)

!----------------------data analysis--beigns-  --------------------------------
  do l=1,N_RAW
     ! 1- sqrt(rho)*u
     f(:,:,:)=dCMPLX(dsqrt(data_in(l,:,:,:,1))*data_in(l,:,:,:,2),0) 
     !f(:,:,:)=dCMPLX(data_in(l,:,:,:,1)*(data_in(l,:,:,:,2)**2+data_in(l,:,:,:,5)**2+data_in(l,:,:,:,6)**2)/2.0d0,0) 
         ! Kinetic energy = rho*0.5*uiui (:,:,:,j) get for each j-time

!%%%%%%%%%%%%%%%%%%%%%%%% test
    !do i=1,N_MESH;do j=1,N_MESH;do k=1,N_MESH;
    !   f(i,j,k)=dcos(2d0*pi*(i)/N_MESH)
    !enddo;enddo;enddo
!%%%%%%%%%%%%%%%%%%%%%%%% test 

    ! -- for every j,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(:,j,k));
        f(:,j,k)=f(:,j,k)/N_MESH
       CALL fftshift(f(:,j,k),N_MESH)
    !write(*,*) "good"
    end do;end do
    ! --

    !-- for every i1,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,:,k));
       f(j,:,k)=f(j,:,k)/N_MESH
       CALL fftshift(f(j,:,k),N_MESH)
    end do;end do
    ! --

    !-- for every i1,i2
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,k,:));
       f(j,k,:)=f(j,k,:)/N_MESH
       CALL fftshift(f(j,k,:),N_MESH)
    end do;end do
  
    fu=f;

    ! 2- sqrt(rho)*v
     f(:,:,:)=dCMPLX(dsqrt(data_in(l,:,:,:,1))*data_in(l,:,:,:,5),0) 

    ! -- for every j,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(:,j,k));
        f(:,j,k)=f(:,j,k)/N_MESH
       CALL fftshift(f(:,j,k),N_MESH)
    !write(*,*) "good"
    end do;end do
    ! --

    !-- for every i1,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,:,k));
       f(j,:,k)=f(j,:,k)/N_MESH
       CALL fftshift(f(j,:,k),N_MESH)
    end do;end do
    ! --

    !-- for every i1,i2
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,k,:));
       f(j,k,:)=f(j,k,:)/N_MESH
       CALL fftshift(f(j,k,:),N_MESH)
    end do;end do

   fv=f;

 ! 3- sqrt(rho)*w
     f(:,:,:)=dCMPLX(dsqrt(data_in(l,:,:,:,1))*data_in(l,:,:,:,6),0) 

    ! -- for every j,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(:,j,k));
        f(:,j,k)=f(:,j,k)/N_MESH
       CALL fftshift(f(:,j,k),N_MESH)
    !write(*,*) "good"
    end do;end do
    ! --

    !-- for every i1,k
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,:,k));
       f(j,:,k)=f(j,:,k)/N_MESH
       CALL fftshift(f(j,:,k),N_MESH)
    end do;end do
    ! --

    !-- for every i1,i2
    do j=1,N_MESH; do k=1,N_MESH
       CALL fft(f(j,k,:));
       f(j,k,:)=f(j,k,:)/N_MESH
       CALL fftshift(f(j,k,:),N_MESH)
    end do;end do

   fw=f;

    ! --

    ! --
    p1=1
    do i=1,N_MESH; do j=1,N_MESH; DO k=1,N_MESH;
       final_index(p1,l)=nint(dsqrt(dble(i-(N_MESH/2+1))**2+ &
                      dble(j-(N_MESH/2+1))**2+dble(k-(N_MESH/2+1))**2))
       final(p1,l)=(abs(fu(i,j,k))**2+abs(fv(i,j,k))**2+abs(fw(i,j,k))**2)/2d0 !*4d0*pi*( & 
                   !dble(i-(N_MESH/2+1))**2 + dble(j-(N_MESH/2+1))**2 + dble(k-(N_MESH/2+1))**2   )/2
       !final(p1,l)=final(p1,l)/N_MESH**3  
       p1=p1+1  
    end do;  end do;  end do
  end do
! --------------- here finish the l=1,N_RAW

  do l=1,n_raw
     nmax(l)=maxval(final_index(:,l))
     nmin(l)=minval(final_index(:,l))
  end do

  nmaxmax=maxval(nmax)
  nminmin=minval(nmin)

  write(*,*) n_mesh,'number'
  allocate(output(nminmin:nmaxmax,1:N_RAW))
  allocate(cum_output(nminmin:nmaxmax,1:N_RAW))
!  write(*,*) nmaxmax

  do l=1,N_raw;do p1=nminmin,nmaxmax
     output(p1,l)=0d0
  end do;end do

  
  do l=1,N_raw;do p1=1,N_MESH**3
     output(final_index(p1,l),l) = output(final_index(p1,l),l) + final(p1,l);
  end do;enddo

do l=1,N_raw;
   do i=nminmin,nmaxmax
   su_output(l)=su_output(l)+output(i,l);
   cum_output(i,l)=su_output(l);
   enddo;
enddo


do l=1,N_RAW
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
 !    if (i.lt.(N_mesh*sqrt(3d0)/2))then ! sqrt(3) -> since our wavenumber is not in 1 direction but a module wavenumber
     if (i.lt.(N_mesh*sqrt(1d0)/2))then 
     if (output(i,l).gt.0.000000000000d0) then
         write(20,*) i,output(i,l),urms2/k0**5*16d0*dsqrt(2d0/pi)*i**4*exp(-2d0*i**2/k0**2)
     !I**4*exp(-2d0*i**2/16)
     endif;end if
  end do
  close(20)

end do
!---------------------box 2: WRITING DATA END----------------------------------------------------
!  end do

!  write(*,*) "good down down"

END PROGRAM H5_MAIN_FFT
