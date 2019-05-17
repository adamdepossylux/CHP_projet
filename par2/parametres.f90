module parametres
integer,parameter :: Nx=400,Ny=400,NxNy=Nx*Ny,Nt=3
double precision,parameter :: Lx=1.d0,Ly=1.d0,dt=0.01d0,D=1.d0,hx=1.d0/(Nx+1),hy=1.d0/(Ny+1)
end module
