	program comp_tools
	
	implicit none
	
	integer i, j, digit, exp_t, arrays_dim, t_slices, t_index
	integer master, i10, i100, lines
	parameter (t_slices=40, arrays_dim=258)
	character *200		header_line
	character *13   	elapsed_t, bounce_t
	character *17   	filename
	character *3   		tail
	double precision 	t_N(t_slices), t_b_N(t_slices), dt
	double precision	u, v, w, shock_t(1000), shock_r(1000)
	double precision	r(t_slices,arrays_dim), dr(t_slices,arrays_dim), q_p
	double precision	lum, rstmss, rho(t_slices,arrays_dim)
	double precision	Temp, d1, d2, d3, d4, d5, d6, d7, d8, lapse(t_slices,arrays_dim)
	double precision	t_proper(t_slices,arrays_dim)
	double precision	enu_lum(t_slices,arrays_dim), e_nu_av(t_slices,arrays_dim)
	double precision	e_nu_rms(t_slices,arrays_dim), e_nu_flx(t_slices,arrays_dim)
	double precision	anu_lum(t_slices,arrays_dim), a_nu_av(t_slices,arrays_dim)
	double precision	a_nu_rms(t_slices,arrays_dim), a_nu_flx(t_slices,arrays_dim)
	double precision	tnu_lum(t_slices,arrays_dim), t_nu_av(t_slices,arrays_dim)
	double precision	t_nu_rms(t_slices,arrays_dim), t_nu_flx(t_slices,arrays_dim)
	double precision	tbnu_lum(t_slices,arrays_dim), tb_nu_av(t_slices,arrays_dim)
	double precision	tb_nu_rms(t_slices,arrays_dim), tb_nu_flx(t_slices,arrays_dim)
	
c	Initialization

	t_index = 1
	
	
c	Master Loop: Goes through all the files

        do master=0,t_slices
	
      		i100     = master/100
      		i10      = (master-100*i100)/10
      		tail     = char(i100+48)//char(i10+48)//char(master-10*i10-100*i100+48)
		filename = "model0001_00"//tail//".d"
		write(*,*) filename  !; pause
	
		open (unit=1,file=filename,status='old')
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Header analysis

c	Get rid of the first lines
		lines = 9
		if (master > 0) lines = 915
		do i=1, lines
	  		read(1,'(A)') header_line
		enddo

c	Recover elapsed time and time from bounce
		read(1,'(A)') header_line
		! write(*,*) header_line
		elapsed_t = header_line(24:36)
		bounce_t  = header_line(66:80)
c	write(*,*) elapsed_t, bounce_t

c	Convert the strings to real
        
c	mantissas
		t_N(t_index) = 0.d0; t_b_N(t_index) = 0.d0
		t_N(t_index) = iachar(elapsed_t(1:1))-48; t_b_N(t_index) = iachar(bounce_t(1:1))-48
	
		do i = 3, 9
	 	  t_N(t_index)   = t_N(t_index)   + (iachar(elapsed_t(i:i))-48) * 10.d0**(2.d0-i)
	 	  t_b_N(t_index) = t_b_N(t_index) + (iachar(bounce_t(i:i))-48) * 10.d0**(2.d0-i)
c	 write(*,*) digit,t_N(t_index)
		enddo

c	exponents
		exp_t = 10.d0 * (iachar(elapsed_t(12:12))-48) + iachar(elapsed_t(13:13))-48
		if (elapsed_t(11:11) == '-') exp_t = -exp_t
		t_N(t_index) = t_N(t_index) * 10.d0**exp_t
	
		exp_t = 10.d0 * (iachar(bounce_t(12:12))-48) + iachar(bounce_t(13:13))-48
		if (bounce_t(11:11) == '-') exp_t = -exp_t
		t_b_N(t_index) = t_b_N(t_index) * 10.d0**exp_t
	
		write(*,100) t_N(t_index), t_b_N(t_index)
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Recovering the data


c	Get rid of the middle lines
		do i=1,33
	  	  read(1,'(A)') header_line
	 	  ! write(*,*) header_line
		enddo
	
c	Read Data "Model Configuration"
		j = 300
		do while (j /= 2)
	  	  read(1,200)  j, u, v, w, r(t_index,j),dr(t_index,j), q_p, lum, rstmss, 
     &			rho(t_index,j), Temp
c	  	  write(*,300) j, u, v, w, r(t_index,j), dr(t_index,j), q_p, lum, rstmss, 
c     &			rho(t_index,j), Temp
		enddo
		
c	Get rid of the middle lines
		lines = 1350
		if (mod(master,5) /= 0) lines = 549
		do i=1, lines
	  	  read(1,'(A)') header_line
c	 	 write(*,*) header_line
		enddo
		
c	Read Data "Hydrodynamic data - Relativistic parameters"
		j = 300
		do while (j /= 2)
	  	  read(1,400)  j, d1, d2, d3, d4, d5, d6, d7, d8, lapse(t_index,j)
c	  	  write(*,500) j, d1, d2, d3, d4, d5, d6, d7, d8, lapse(t_index,j)

		enddo
	
	if (mod(master,5)==0) then
	
c	Get rid of the middle lines
		lines = 2982
		if (master == 0 )  lines = 3249
		if (master > 10 ) lines = 2993
		do i=1,lines
	  	  read(1,'(A)') header_line
c	 	 write(*,*) header_line
		enddo
		
c	Read Data "electron neutrino energy data"
		j = 300
		do while (j /= 2)
	  	  read(1,600)  j, d1, d2, d3, d4, d5, d6, d7,
     &			enu_lum(t_index,j), e_nu_av(t_index,j), e_nu_rms(t_index,j), e_nu_flx(t_index,j)
c	  	  write(*,700) j, d1, d2, d3, d4, d5, d6, d7,
c     &			enu_lum(t_index,j), e_nu_av(t_index,j), e_nu_rms(t_index,j), e_nu_flx(t_index,j)
		enddo
		
c	Get rid of the middle lines
		do i=1,11
	  	  read(1,'(A)') header_line
c	 	 write(*,*) header_line
		enddo
		
c	Read Data "electron antineutrino energy data"
		j = 300
		do while (j /= 2)
	  	  read(1,600)  j, d1, d2, d3, d4, d5, d6, d7,
     &			anu_lum(t_index,j), a_nu_av(t_index,j), a_nu_rms(t_index,j), a_nu_flx(t_index,j)
c	  	  write(*,700) j, d1, d2, d3, d4, d5, d6, d7,
c     &			anu_lum(t_index,j), a_nu_av(t_index,j), a_nu_rms(t_index,j), a_nu_flx(t_index,j)
		enddo
		
c	Get rid of the middle lines
		do i=1,11
	  	  read(1,'(A)') header_line
c	 	 write(*,*) header_line
		enddo
		
c	Read Data "x-neutrino energy data"
		j = 300
		do while (j /= 2)
	  	  read(1,600)  j, d1, d2, d3, d4, d5, d6, d7,
     &			tnu_lum(t_index,j), t_nu_av(t_index,j), t_nu_rms(t_index,j), t_nu_flx(t_index,j)
c	  	  write(*,700) j, d1, d2, d3, d4, d5, d6, d7,
c     &			tnu_lum(t_index,j), t_nu_av(t_index,j), t_nu_rms(t_index,j), t_nu_flx(t_index,j)
		enddo
		
c	Get rid of the middle lines
		do i=1,11
	  	  read(1,'(A)') header_line
c	 	 write(*,*) header_line
		enddo
		
c	Read Data "x-antineutrino energy data"
		j = 300
		do while (j /= 2)
	  	  read(1,600)  j, d1, d2, d3, d4, d5, d6, d7,
     &			tbnu_lum(t_index,j), tb_nu_av(t_index,j), tb_nu_rms(t_index,j), tb_nu_flx(t_index,j)
c	  	  write(*,700) j, d1, d2, d3, d4, d5, d6, d7,
c     &			tbnu_lum(t_index,j), tb_nu_av(t_index,j), tb_nu_rms(t_index,j), tb_nu_flx(t_index,j)
		enddo
	
	endif  ! End of if (mod(master,5)==0)
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	t_index = t_index + 1
	
	enddo    ! Master loop	
	
	close(unit=1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Read shock position

	open (unit=1,file='shock.d',status='old')
	
	do j = 1, 337
	  read(1,800)  shock_t(j), shock_r(j)
	  ! write(*,900) shock_t(j), shock_r(j)
	enddo
		
	close(unit=1)

	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Calculate proper time at the origin

	print *, "  "
	print *, "  Rest-mass Density at the center"
        print *, "	time	     Proper time	rho"
	print *, "  "
	
	do master = 0, t_slices

	  if (master == 0) then
	    t_proper(master+1,2) = lapse(master+1,2) * t_N(master+1)
	  else
            dt = t_N(master+1) - t_N(master)
	    t_proper(master+1,2) = lapse(master+1,2) * dt + t_proper(master,2)
	  endif
	  	
           
c	  write(*,33) master,r(master+1,2),lapse(master+1,2)
	  write(*,34) t_N(master+1),t_proper(master+1,2),rho(master+1,2)
	  
	enddo ! master
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	Calculate proper time at the shock

	print *, "  "
	print *, "  Rest-mass Density at the shock"
        print *, "	time	     Proper time	rho"
	print *, "  "
	
	do master = 0, t_slices
          
	  do j = 2, 257
	        write(*,*) master, j, r(master+1,j), shock_r(j)
	  	if (r(master+1,j) == shock_r(j)) goto 22
	  enddo

22	  continue

	  if (master == 0) then
	    t_proper(master+1,2) = lapse(master+1,2) * t_N(master+1)
	  else
            dt = t_N(master+1) - t_N(master)
	    t_proper(master+1,2) = lapse(master+1,2) * dt + t_proper(master,2)
	  endif
	  	
           
c	  write(*,33) master,r(master+1,2),lapse(master+1,2)
	  write(*,34) t_N(master+1),t_proper(master+1,2),rho(master+1,2)
	  
	enddo ! master
	
33	format(i6,f14.0,f14.4)
34	format(2f14.8,e16.4)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
100	format (3e16.8)
200	format (i5,6e11.3, e12.3, e14.6, 2e11.3)
300	format (i5,6e12.4, e13.4, e15.7, 2e12.4)
400	format (i4,9e11.3)
500	format (i5,9e12.4)
600	format (i4,11e11.4)
700	format (i5,11e12.5)
800	format (e15.8,e12.4)
900	format (e16.9,e13.5)

	end program comp_tools
