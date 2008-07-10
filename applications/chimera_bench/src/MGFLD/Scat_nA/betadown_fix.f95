ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      function beta_down_mgfld
c
c
c      Purpose: to compute the kernal for down transitions in inelastic neutrino-nucleus
c                   interactions, as discussed in Fuller and Meyer 1991, ApJ, 376, 701-716.
c                   These kernels are for use in exothermic scattering (eq. 2.7 of
c                   Fuller and Meyer) and nuclear de-excitation into pairs (eq. 2.8a).
c                   This routine is for MGFLD.  The angle dependence is integrated out.
c      Input: 
c               z = atomic number of the nucleus
c               a = mass number of the nucleus
c               tmev = matter temperature in MeVs
c               e_representative = representative nucleus energy in the thermal bath, called
c                                  <E> in Fuller and Meyer.  <E> = aT^2, where T is the
c                                  local matter temperature.  a is typically given by
c                                  a = A/8 MeV^-1
c               delta_e = change in the energy of the nucleus during the transition.  Note
c                            delta_e = e_neutrino_out - e_neutrino_in for exothermic
c                            scattering and delta_e = e_neutrino_out + e_neutrino_bar_out
c                            for de-excitation into pairs
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function beta_down_mgfld( z, a, tmev,
     &         e_representative, delta_e)
c
      implicit none
      double precision z, a, tmev, e_representative, delta_e, mu_p,
     &       mu_n
      double precision f1, f2, e_1, e_2, e_3, e_4
      double precision lambda_1, lambda_2, lambda_3, lambda_4
      double precision f
      double precision p1, p2, p3, I_0_down
      double precision q1, q2, I_A
      double precision s0, r0, s1, s2, hbarc, I_F

c     Fit parameters (Table 2 of Fuller and Meyer 1991)

      data f1, f2, e_1, lambda_1, e_2, lambda_2, e_3, lambda_3, e_4,
     &        lambda_4/0.09d0, 0.017d0, 3.0d0, 1.0d0, 4.0d0, 1.0d0,
     &        17.0d0, 8.0d0, 22.0d0, 1.0d0/

c     Nuclear parameters

      data r0, hbarc/1.2d0, 197.d0/

c     Compute chemical potentials

      mu_p = 50.d0*(z/30.d0)**(2.d0/3.d0)
      mu_n = 50.d0*((a-z)/30.d0)**(2.d0/3.d0)

c     Compute f

      f = f1 + f2*tmev**2

c     Compute I_0_down

      if(delta_e.ne.0.d0) then
         p1 = 79.d0*f/(dexp(delta_e/tmev) -1.d0)
         p2 = (z/30.d0)*((1.d0+delta_e/mu_p)**(3.d0/2.d0) - 1.d0)
         p3 = ((a-z)/30.d0)*((1.d0+delta_e/mu_n)**(3.d0/2.d0) - 1.d0)
         I_0_down = p1*(p2 + p3)
      else
         I_0_down = 79d0*f*1.5d0*tmev*((z/30.d0)/mu_p +
     &                  ((a-z)/30.d0)/mu_n)
      endif

c     Compute I_A

      q1 = 1.d0/(dexp((e_1 - delta_e)/lambda_1) + 1.d0)
      q2 = 1.d0/(dexp((delta_e - e_2)/lambda_2) + 1.d0)
      I_A = q1*q2

c     Compute I_F

      s0 = (r0*a**(1.d0/3.d0)*delta_e/hbarc)**2           
      s1 = 1.d0/(dexp((e_3 - delta_e)/lambda_3) + 1.d0)
      s2 = 1.d0/(dexp((delta_e - e_4)/lambda_4) + 1.d0)
      I_F = s0*s1*s2

c     Compute beta_down

      beta_down_mgfld = I_0_down * (I_A + I_F)

      return
      end
