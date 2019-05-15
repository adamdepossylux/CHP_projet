		module fonctions
			implicit none
			contains

				function f1(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: f1
					f1=2*(y-y**2+x-x**2)
				end function f1

				function f(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: f
					f=sin(x)+cos(y)
				end function f

				function g(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: g
					g=0.d0
				end function g

				function h(x,y)
					implicit none
					double precision,intent(in)::x,y
					double precision :: h
					h=1.d0
				end function h




	! fonction f dépendante du temps

				function ff(x,y,t,Ly,Lx)
					implicit none
					double precision, parameter :: pi=3.14159265
					double precision, intent(in)::x,y,t,Ly,Lx
					double precision::ff
					ff=exp(-((x-Lx/2)**2))*(exp(-((y-Ly/2)**2)))*cos((pi/2)*t)
				end function ff


	!Solution exacte


		end module
