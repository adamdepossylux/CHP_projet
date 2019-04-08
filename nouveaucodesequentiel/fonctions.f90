		module fonctions
			implicit none
			contains
				function f(x,y) 
					implicit none
					real, intent(in)::x,y
					real :: f
					f=2*(y-y**2 + x + x**2)
				end function f 
				
				function g(x,y) 
					implicit none
					real, intent(in)::x,y 
					real :: g
					g=0
				end function g
			
				function h(x,y) 
					implicit none
					real,intent(in)::x,y 
					real :: h
					h=0
				end function h
				
				function ff(x,y,t,Ly,Lx) 
					implicit none
					real, parameter :: pi=3.14159265
					real, intent(in)::x,y,t,Ly,Lx
					real::ff
					ff=exp(-((x-Lx/2)**2))*(exp(-((y-Ly/2)**2)))*cos((pi/2)*t)
				end function ff
				
				function gg(x,y,t)
					implicit none
					real, intent(in)::x,y,t
					real :: gg
					gg=0
				end function gg
			
				
				function hh(x,y,t) 
					implicit none
					real, intent(in)::x,y,t
					real :: hh
					hh=1.
				end function hh

				
		end module
			