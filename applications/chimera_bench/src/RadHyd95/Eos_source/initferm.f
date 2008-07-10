C23456789012345678901234567890123456789012345678901234567890123456789012
C
c      Program to compute spline fits to fermi integrals
cc  Must provide data file 94
      subroutine initferm
      implicit  double precision (a-h,o-z)
      parameter (n=201)
      dimension f32(n),f12(n),fm12(n),eta(n),fr(n)
      dimension f32a(n),f12a(n),fra(n),fia(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
C
      open(94,file='../../Equation_of_State/Eos_lsu/fermi.atb',
     & status='old')
C
      do 10 i=1,n
       read(94,*)eta(i),f32(i),f12(i),fm12(i)
 10    fr(i)=f12(i)/fm12(i)
C
      close(94,status='keep')
C
       call spline(eta,f12,n,f12a)
c       write(*,1)f12a
       call spline(eta,f32,n,f32a)
c       write(*,1)f32a
       call spline(eta,fr,n,fra)
c       write(*,1)fra
       call spline(f12,eta,n,fia)
c       write(*,1)fia
 1     format(8(1pe10.3))
       return
       end
