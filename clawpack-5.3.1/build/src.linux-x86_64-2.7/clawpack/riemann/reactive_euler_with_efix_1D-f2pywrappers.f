C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pyinitcparam(setupfunc)
      external setupfunc
      double precision gamma
      double precision qheat
      common /cparam/ gamma,qheat
      call setupfunc(gamma,qheat)
      end

