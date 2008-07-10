      subroutine eos0
c***********************************************************************
      USE array_module
      USE eos_bck_module
      USE physcnst_module
      USE edit_module, ONLY : nlog
      implicit double precision (a-h,o-z)
      logical first
      save
c***********************************************************************
      parameter(third=1./3.,ba=-7.075,zero=0.0,dnp=1.2935)
      parameter(a0=.067,c0=2.3798e-4)
c                                                                      c
      data first/.true./
c***********************************************************************
      f3(z)        = sign( (abs(z))**third , z)
      fw1(x)       = 78.0 * x * x * (1.0-x) * f3(1.0-x)
      fbulk(x)     = wnm + ws * (1.0-2.0*x)*(1.0-2.0*x)
      fphi(x)      = 1.0 - 3.0 * (0.5-x)**2
      fxk(zabck)   = xk0*(1. - xkzafac*(1.- 2.*zabck)**2)
      fcomp(theta) = (fxk(zabck)/18.0) * (1.0 - theta)*(1.0 - theta)
c******************************************** get energy zeroes ********
              if(first)then
      d00 = 0.16
      yeFe = 26./56.
      zabck = yeFe
      theta0 = 1.0
      size0 = fw1(yeFe) /(fphi(yeFe))**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      theta0 = 1 +(3.0/fxk(yeFe))*size0/theta0/theta0**third
      size0 = size0/theta0**third
      egy0 = fbulk(yeFe) + size0 + fcomp(theta0)
      write(nlog,999)yeFe,d00*theta0*fphi(yeFe),theta0,size0,egy0
      first = .false.
             endif
c***********************************************************************
                        return
c***********************************************************************
999   format(' zero level for fe: yebck,d0,theta,wsize,egy',5f10.5)
                        end
