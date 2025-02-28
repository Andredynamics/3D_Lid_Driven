function [Nx,Ny,Nz,x,y,z,h,hquadro,Lx,Ly,Lz,omega] = Params()
Nx = 32; Ny = 32; Nz = 32; 
Lx = 1;   Ly = 1;   Lz = 1;
omega = 1;
z = linspace(0,Lz,Nz);
y = linspace(0,Ly,Ny);
x = linspace(0,Lx,Nx);
h = x(2) - x(1); hquadro = h^2; 

end

