clc; clear; close all; 
[Nx,Ny,Nz,x,y,z,h,hquadro] = Params();
[U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
 V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
 W_Up,W_Down,W_Front,W_Back,W_Right,W_Left,Uref]= BCs();
global Re

Re_vec = [100,400,1000,3200]; it = 1;
C_vec = [0.1,0.1,0.05,0.05];
for cont = 1 : 4 
T = 60; Re = Re_vec(cont); C = C_vec(cont); beta = 0.5; 
Uref = 1;
dt = min([C*h/Uref,beta*hquadro*Re]); 
Nt = round(T/dt);
dt = T / Nt;
t = linspace(0,T,Nt);

U = zeros(Ny - 1, Nx, Nz - 1); 
V = zeros(Ny, Nx - 1, Nz - 1); 
W = zeros(Ny - 1, Nx - 1, Nz); 
%% Definisco numgrid, laplaciano e BCs su laplaciano. 
G = zeros(Ny + 1, Nx + 1, Nz + 1); 
num = 0;
for k = 2 : Nz 
    for j = 2 : Nx
        for i = 2 : Ny
            num = num + 1; G(i,j,k) = num;
        end
    end
end

Lap = Delsq_3D_Lap(G) / hquadro; 

for i = 2 : Ny 
    for j = 2 : Nx
        k = 2;  num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro; 
        k = Nz; num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro;
    end
end

for i = 2 : Ny 
    for k = 2 : Nz
        j = 2;  num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro; 
        j = Nx; num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro; 
    end
end

for j = 2 : Nx
    for k = 2 : Nz 
        i = 2;  num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro; 
        i = Ny; num = G(i,j,k); Lap(num,num) = Lap(num,num) + 1 / hquadro; 
    end
end
%% Comincio iterazione. 
tol = 1e-5;
maxvar = 1;
while maxvar > tol 
    U1 = U; V1 = V; W1 = W; [FU1,FV1,FW1] = RHS(U1,V1,W1);

    U2 = U + 0.5*dt*FU1; V2 = V + 0.5*dt*FV1; W2 = W + 0.5*dt*FW1;
    Div = Div_Calc(U2,V2,W2); Div = Div(:);
    P = Lap\Div; P = reshape(P,Ny-1,Nx-1,Nz-1);
    [Px,Py,Pz] = Grad_Calc(P);
    U2 = U2 - Px; V2 = V2 - Py; W2 = W2 - Pz; 
    [FU2,FV2,FW2] = RHS(U2,V2,W2);

    U3 = U + 0.5*dt*FU2; V3 = V + 0.5*dt*FV2; W3 = W + 0.5*dt*FW2;
    Div = Div_Calc(U3,V3,W3); Div = Div(:);
    P = Lap\Div; P = reshape(P,Ny-1,Nx-1,Nz-1);
    [Px,Py,Pz] = Grad_Calc(P);
    U3 = U3 - Px; V3 = V3 - Py; W3 = W3 - Pz; 
    [FU3,FV3,FW3] = RHS(U3,V3,W3);

    U4 = U + dt*FU3; V4 = V + dt*FV3; W4 = W + dt*FW3;
    Div = Div_Calc(U4,V4,W4); Div = Div(:);
    P = Lap\Div; P = reshape(P,Ny-1,Nx-1,Nz-1);
    [Px,Py,Pz] = Grad_Calc(P);
    U4 = U4 - Px; V4 = V4 - Py; W4 = W4 - Pz; 
    [FU4,FV4,FW4] = RHS(U4,V4,W4);
    
    U = U + ((1/6)*(FU1+FU4)+(1/3)*(FU2+FU3))*dt;
    V = V + ((1/6)*(FV1+FV4)+(1/3)*(FV2+FV3))*dt;
    W = W + ((1/6)*(FW1+FW4)+(1/3)*(FW2+FW3))*dt;
    Div = Div_Calc(U,V,W); Div = Div(:);
    P = Lap\Div; P = reshape(P,Ny-1,Nx-1,Nz-1);
    [Px,Py,Pz] = Grad_Calc(P);
    U = U - Px; V = V - Py; W = W - Pz; 
    
    %disp(num2str(max(max(max(abs(Div_Calc(U,V,W)))))));
%% Plotto la divergenza max e la variazione della soluzione.     
    [X,Y,Z] = meshgrid(x(1:Nx-1),y(1:Ny-1),z(1:Nz-1));
    
    if it ~= 1
    Iu_cen_old = Iu_cen; 
    Iv_cen_old = Iv_cen;
    Iw_cen_old = Iw_cen; 
    end

    Iu_cen = zeros(Ny - 1, Nx - 1, Nz - 1);
    Iv_cen = zeros(Ny - 1, Nx - 1, Nz - 1); 
    Iw_cen = zeros(Ny - 1, Nx - 1, Nz - 1);
    for i = 1 : Ny - 1
       for j = 1 : Nx - 1
          for k = 1 : Nz - 1
             Iu_cen(i,j,k) = (U(i,j + 1,k) + U(i,j,k))/2;
             Iv_cen(i,j,k) = (V(i + 1,j,k) + V(i,j,k))/2;
             Iw_cen(i,j,k) = (W(i,j,k + 1) + W(i,j,k))/2;
          end
       end
    end
   if it ~= 1 
      umaxvar = max(max(max(abs(Iu_cen_old - Iu_cen))));
      vmaxvar = max(max(max(abs(Iv_cen_old - Iv_cen))));
      wmaxvar = max(max(max(abs(Iw_cen_old - Iw_cen))));
      maxvar = max([umaxvar,vmaxvar,wmaxvar])/(dt*Re);
      disp(num2str(maxvar));
   end
   
   if it ~= 1
   time = horzcat(time,t(it)); %#ok<AGROW>
   Div_plot = horzcat(Div_plot,max(max(max(abs(Div_Calc(U,V,W)))))); %#ok<AGROW>
   Var_plot = horzcat(Var_plot,maxvar); %#ok<AGROW>
   else
   time = t(it);
   Div_plot = max(max(max(abs(Div_Calc(U,V,W))))); 
   Var_plot = 0;
   end
   

   it = it + 1; 

   figure(2);
   plot(time,Div_plot,'.-k'); hold on; drawnow; 
   
   figure(3);
   if it ~= 1
   plot(time,Var_plot,'.-k'); hold on; drawnow; 
   end

   % figure(3);
   % title('Streamslice lungo z');
   % for num = 1 : 8
   % subplot(3,3,num); axis square;
   % U_2D = squeeze(Iu_cen(:,:,round(num/9*Nz)));
   % V_2D = squeeze(Iv_cen(:,:,round(num/9*Nz)));
   % Psi = UVW2Psi(U_2D,V_2D);
   % contour(y(2:Nx-1),x(2:Nx-1),Psi);
   % end
   % 
   % figure(4);
   % title('Streamslice lungo y');
   % for num = 1 : 8 
   % subplot(3,3,num); axis square;
   % U_2D = squeeze(Iu_cen(round(num/9*Ny),:,:));
   % W_2D = squeeze(Iw_cen(round(num/9*Ny),:,:));
   % Psi = UVW2Psi(U_2D,W_2D);
   % contour(y(2:Nx-1),x(2:Nx-1),Psi);
   % end
   % 
   % figure(5);
   % title('Streamslice lungo x');
   % for num = 1 : 8
   % subplot(3,3,num); axis square;
   % V_2D = squeeze(Iu_cen(:,round(num/9*Nx),:));
   % W_2D = squeeze(Iw_cen(:,round(num/9*Nx),:));
   % Psi = UVW2Psi(V_2D,W_2D);
   % contour(y(2:Nx-1),z(2:Nx-1),Psi);
   % end
end

close all;
end
