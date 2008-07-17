! find triangle which contaion given point

subroutine find_triangle(r,z,it,p)
  implicit none
  real (kind=8) :: r,z,p(3)
  integer :: i,j,it
  integer, external :: is_inside

    
  !check vaildity of range
  if( r>= r_max .or. r<=r_min .or. z>=z_max .or. z<=z_min) then
     print *, 'error in find_triangle',r,z
     print *, r_min, r_max, z_min, z_max
     stop
  endif

  ! find index in a rectanglar grid  
  i= (r- r_min)/NR
  j= (z- z_min)/NZ
  
  ! search triangular
  find=0;k=1
  do while(find/=1)
     !error eheck
     if(k > n_triangle(i,j,0) ) then
        print *, 'error in find_triangle',r,z
        print *, 'no triangle found'
        stop
     endif
     
     !get triangle number
     call get_coeff( r,z,n_triangle(i,j,k),p)
     if( is_inside(p) ) then
        find=1
        it=n_triangle(i,j,k)
     endif
     k=k+1
  enddo

end subroutine find_triangle

subroutine init_search
  implicit none
  
  ! for all rectangle
  do i=1, NR
     do j=1, NZ
        n=0
        ! for all triangle
        do k=1, NT
           overlapping=0
           
           ! Node of T inside R test
           !for all 3 nodes of T
           do nd=1, 3
              !check inside R
              r=node( triangle_node(k,nd) )%r
              z=node( triangle_node(k,nd) )%z
              !if inside, overlapping=1 and exit
              if( is_inside_R(r,z,i,j)==1 ) then
                 overlapping=1
                 exit
              endif
           enddo

           ! T edge - R edge intersection test
           if(overlapping/=1) then
              do edge_set=1, 12
                 tedge= edge_set/4
                 redge= mode(edge_set,4)+1
                 
                 ! check intersection
                 

                    
                 ! if interect, overlapping=1 and exit
                 if( is_intersect(l1,l2) ) then
                    overlapping=1
                    exit
                 endif
              enddo
           endif
           
           ! triangle inside rectangle test
           if(overlapping/=1) then
              


           endif
           !if it is found, store the triangle number and n=n+1
        enddo
        n_triangle(i,j,0)=n
     enddo
  enddo


end subroutine initilizing
