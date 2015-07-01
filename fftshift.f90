subroutine fftshift(f,N)
    USE fft_mod
    implicit none
    integer :: N
    COMPLEX(kind=dp):: f(N),ftemp(N)
    
    ftemp=f
    f(1:N/2)=ftemp(N/2+1:N)
    f(N/2+1:N)=ftemp(1:N/2)

    
return
end 
    
