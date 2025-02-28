function [U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
          V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
          W_Up,W_Down,W_Front,W_Back,W_Right,W_Left,Uref]= BCs()

[Nx,Ny,Nz,x,y,z,h,hquadro,Lx,Ly,Lz,omega] = Params();
% 
% [X, Z] = meshgrid(x, z);

% raggi = sqrt((X - Lx/2).^2 + (Z - Lz/2).^2);
% angoli = atan2(Z - Lz/2, X - Lx/2);
% vtan = omega * (2 * pi * raggi);
% u = vtan .* cos(angoli - pi/2);  
% w = vtan .* sin(angoli - pi/2);  
% quiver(X, Z, u, w);

Uref = 1;
U_Up = 1; 
U_Down = 0; 
U_Front = 0; 
U_Back = 0; 
U_Right = 0; 
U_Left = 0; 

V_Up = 0;
V_Down = 0; 
V_Front = 0; 
V_Back = 0; 
V_Right = 0; 
V_Left = 0; 

W_Up = 0; 
W_Down = 0; 
W_Front = 0; 
W_Back = 0; 
W_Right = 0; 
W_Left = 0; 
end

