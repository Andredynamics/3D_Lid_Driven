function [FU,FV,FW] = RHS(U,V,W)
[Nx,Ny,Nz,~,~,~,h,hquadro,~,~,~,~] = Params();
global Re
[U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
 V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
 W_Up,W_Down,W_Front,W_Back,W_Right,W_Left]= BCs();

FU = zeros(Ny - 1, Nx, Nz - 1);
FV = zeros(Ny, Nx - 1, Nz - 1);
FW = zeros(Ny - 1, Nx - 1, Nz);

%% Interpolo sui centri e sui centri dei lati. 
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

Iu_xz = zeros(Ny - 1, Nx, Nz);
Iw_xz = zeros(Ny - 1, Nx, Nz);
Iu_xy = zeros(Ny, Nx, Nz - 1);
Iv_xy = zeros(Ny, Nx, Nz - 1);
Iv_yz = zeros(Ny, Nx - 1, Nz);
Iw_yz = zeros(Ny, Nx - 1, Nz);

i = 1 : Ny - 1; j = 2 : Nx - 1; k = 2 : Nz - 1; 
Iu_xz(i,j,k) = (U(i,j,k) + U(i,j,k - 1))/2;
Iw_xz(i,j,k) = (W(i,j,k) + W(i,j - 1,k))/2;

i = 2 : Ny - 1; j = 1 : Nx - 1; k = 2 : Nz - 1; 
Iv_yz(i,j,k) = (V(i,j,k) + V(i,j,k - 1))/2;
Iw_yz(i,j,k) = (W(i,j,k) + W(i - 1,j,k))/2;

i = 2 : Ny - 1; j = 2 : Nx - 1; k = 1 : Nz - 1;
Iu_xy(i,j,k) = (U(i,j,k) + U(i - 1,j,k))/2;
Iv_xy(i,j,k) = (V(i,j,k) + V(i,j - 1,k))/2;

for i = 1 : Ny - 1
    for j = 1 : Nx
        k = 1; Iu_xz(i,j,k) = U_Front; Iw_xz(i,j,k) = W_Front; 
        k = Nz; Iu_xz(i,j,k) = U_Back; Iw_xz(i,j,k) = W_Back;
    end
end

for i = 1 : Ny - 1
    for k = 1 : Nz
        j = 1; Iu_xz(i,j,k) = U_Left; Iw_xz(i,j,k) = W_Left; 
        j = Nx; Iu_xz(i,j,k) = U_Right; Iw_xz(i,j,k) = W_Right;
    end
end

for j = 1 : Nx - 1
     for i = 1 : Ny
         k = 1; Iv_yz(i,j,k) = V_Front; Iw_yz(i,j,k) = W_Front;
         k = Nz; Iv_yz(i,j,k) = V_Back; Iw_yz(i,j,k) = W_Back;
     end
end

for j = 1 : Nx - 1
    for k = 1 : Nz
        i = 1; Iv_yz(i,j,k) = V_Up; Iw_yz(i,j,k) = W_Up;
        i = Ny; Iv_yz(i,j,k) = V_Down; Iw_yz(i,j,k) = W_Down; 
    end
end

for k = 1 : Nz - 1
    for j = 1 : Nx 
        i = 1; Iu_xy(i,j,k) = U_Up; Iv_xy(i,j,k) = V_Up; 
        i = Ny; Iu_xy(i,j,k) = U_Down; Iv_xy(i,j,k) = V_Down;
    end
end


for k = 1 : Nz - 1
    for i = 1 : Ny
        j = 1; Iu_xy(i,j,k) = U_Left; Iv_xy(i,j,k) = V_Left;
        j = Nx; Iu_xy(i,j,k) = U_Right; Iv_xy(i,j,k) = V_Right; 
    end
end

Iu_cen_q = Iu_cen.*Iu_cen; 
Iv_cen_q = Iv_cen.*Iv_cen;
Iw_cen_q = Iw_cen.*Iw_cen;

Iuv = Iu_xy.*Iv_xy;
Iuw = Iu_xz.*Iw_xz;
Ivw = Iv_yz.*Iw_yz;

%% Calcolo termine convettivo. 
duudx = zeros(Ny - 1, Nx, Nz - 1); 
dvvdy = zeros(Ny, Nx - 1, Nz - 1); 
dwwdz = zeros(Ny - 1, Nx - 1, Nz);

i = 1 : Ny - 1; j = 2 : Nx - 1; k = 1 : Nz - 1;
duudx(i,j,k) = (Iu_cen_q(i,j,k) - Iu_cen_q(i,j - 1,k))/h;

i = 2 : Ny - 1; j = 1 : Nx - 1; k = 1 : Nz - 1; 
dvvdy(i,j,k) = (Iv_cen_q(i,j,k) - Iv_cen_q(i - 1,j,k))/h;

i = 1 : Ny - 1; j = 1 : Nx - 1; k = 2 : Nz - 1; 
dwwdz(i,j,k) = (Iw_cen_q(i,j,k) - Iw_cen_q(i,j,k - 1))/h;

duvdy = zeros(Ny - 1, Nx, Nz - 1); 
duwdz = zeros(Ny - 1, Nx, Nz - 1);
i = 1 : Ny - 1; j = 2 : Nx - 1; k = 1 : Nz - 1; 
duvdy(i,j,k) = (Iuv(i + 1,j,k) - Iuv(i,j,k))/h;
duwdz(i,j,k) = (Iuw(i,j,k + 1) - Iuw(i,j,k))/h;

duvdx = zeros(Ny, Nx - 1, Nz - 1); 
dvwdz = zeros(Ny, Nx - 1, Nz - 1);
i = 2 : Ny - 1; j = 1 : Nx - 1; k = 1 : Nz - 1; 
duvdx(i,j,k) = (Iuv(i,j + 1,k) - Iuv(i,j,k))/h;
dvwdz(i,j,k) = (Ivw(i,j,k + 1) - Ivw(i,j,k))/h;

dwudx = zeros(Ny - 1, Nx - 1, Nz); 
dwvdy = zeros(Ny - 1, Nx - 1, Nz);
i = 1 : Ny - 1; j = 1 : Nx - 1; k = 2 : Nz - 1;
dwudx(i,j,k) = (Iuw(i,j + 1,k) - Iuw(i,j,k))/h;
dwvdy(i,j,k) = (Ivw(i + 1,j,k) - Ivw(i,j,k))/h;

FU = FU - duudx - duvdy - duwdz;
FV = FV - dvvdy - duvdx - dvwdz; 
FW = FW - dwwdz - dwudx - dwvdy; 

%% Calcolo termine viscoso. 

LapU = zeros(Ny - 1, Nx, Nz - 1); 
LapV = zeros(Ny, Nx - 1, Nz - 1); 
LapW = zeros(Ny - 1, Nx - 1, Nz); 

i = 2 : Ny - 2; j = 2 : Nx - 1; k = 2 : Nz - 2;
LapU(i,j,k) = (1/Re)*(-6*U(i,j,k) + U(i + 1,j,k) + U(i - 1,j,k) + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + U(i,j,k - 1))/hquadro;

i = 2 : Ny - 1; j = 2 : Nx - 2; k = 2 : Nz - 2;
LapV(i,j,k) = (1/Re)*(-6*V(i,j,k) + V(i + 1,j,k) + V(i - 1,j,k) + ...
                      + V(i,j + 1,k) + V(i,j - 1,k) + V(i,j,k + 1) + ...
                      + V(i,j,k - 1))/hquadro;

i = 2 : Ny - 2; j = 2 : Nx - 2; k = 2 : Nz - 1;
LapW(i,j,k) = (1/Re)*(-6*W(i,j,k) + W(i + 1,j,k) + W(i - 1,j,k) + ...
                      + W(i,j + 1,k) + W(i,j - 1,k) + W(i,j,k + 1) + ...
                      + W(i,j,k - 1))/hquadro;

% LapU: 
for j = 2 : Nx - 1
    for i = 2 : Ny - 2
        k = 1; 
        LapU(i,j,k) = (1/Re)*(-7*U(i,j,k) + U(i + 1,j,k) + U(i - 1,j,k) + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + 2*U_Front)/hquadro;
        k = Nz - 1;
        LapU(i,j,k) = (1/Re)*(-7*U(i,j,k) + U(i + 1,j,k) + U(i - 1,j,k) + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k - 1) + ...
                      + 2*U_Back)/hquadro;
    end
end

for j = 2 : Nx - 1
    for k = 2  : Nz - 2
        i = 1; 
        LapU(i,j,k) = (1/Re)*(-7*U(i,j,k) + U(i + 1,j,k) + 2*U_Up + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + U(i,j,k - 1))/hquadro;
        i = Ny - 1;
        LapU(i,j,k) = (1/Re)*(-7*U(i,j,k) + U(i - 1,j,k) + 2*U_Down + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + 2*U(i,j,k - 1))/hquadro;
    end
end

for j = 2 : Nx - 1
    i = 1; k = 1; 
        LapU(i,j,k) = (1/Re)*(-8*U(i,j,k) + U(i + 1,j,k) + 2*U_Up + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + 2*U_Front)/hquadro;
    i = Ny - 1; k = 1;
        LapU(i,j,k) = (1/Re)*(-8*U(i,j,k) + U(i - 1,j,k) + 2*U_Down + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k + 1) + ...
                      + 2*U_Front)/hquadro;
    i = Ny - 1; k = Nz - 1; 
        LapU(i,j,k) = (1/Re)*(-8*U(i,j,k) + U(i - 1,j,k) + 2*U_Down + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k - 1) + ...
                      + 2*U_Back)/hquadro;
    i = 1; k = Nz - 1; 
        LapU(i,j,k) = (1/Re)*(-8*U(i,j,k) + U(i + 1,j,k) + 2*U_Up + ...
                      + U(i,j + 1,k) + U(i,j - 1,k) + U(i,j,k - 1) + ...
                      + 2*U_Back)/hquadro;
end

% LapV: 
for i = 2 : Ny - 1
    for j = 2 : Nx - 2
        k = 1; 
        LapV(i,j,k) = (1/Re)*(-7*V(i,j,k) + V(i + 1,j,k) + V(i - 1, j, k) + ...
                      + V(i, j + 1, k) + V(i, j - 1, k) + V(i, j, k + 1) + ...
                      + 2*V_Front)/hquadro;
        k = Nz - 1; 
        LapV(i,j,k) = (1/Re)*(-7*V(i,j,k) + V(i + 1,j,k) + V(i - 1, j, k) + ...
                      + V(i, j + 1, k) + V(i, j - 1, k) + V(i, j, k - 1) + ...
                      + 2*V_Back)/hquadro;
    end
end

for i = 2 : Ny - 1
    for k = 2 : Nz - 2
        j = 1; 
        LapV(i,j,k) = (1/Re)*(-7*V(i,j,k) + V(i + 1,j,k) + V(i - 1, j, k) + ...
                      + V(i, j + 1, k) + 2*V_Left + V(i, j, k - 1) + ...
                      + V(i,j,k + 1))/hquadro;
        j = Nx - 1; 
        LapV(i,j,k) = (1/Re)*(-7*V(i,j,k) + V(i + 1,j,k) + V(i - 1, j, k) + ...
                      + V(i, j - 1, k) + 2*V_Right + V(i, j, k - 1) + ...
                      + V(i,j,k + 1))/hquadro;
    end
end

for i = 2 : Ny - 1
    j = 1; k = 1;
        LapV(i,j,k) = (1/Re)*(-8*V(i,j,k) + V(i + 1,j,k) + V(i - 1,j,k) + ...
                      + V(i,j + 1,k) + 2*V_Left + V(i,j, k + 1) + ...
                      + 2*V_Front)/hquadro;
    j = Nx - 1; k = 1; 
        LapV(i,j,k) = (1/Re)*(-8*V(i,j,k) + V(i + 1,j,k) + V(i - 1,j,k) + ...
                      + V(i,j - 1,k) + 2*V_Right + V(i,j, k + 1) + ...
                      + 2*V_Front)/hquadro;
    j = Nx - 1; k = Nz - 1; 
        LapV(i,j,k) = (1/Re)*(-8*V(i,j,k) + V(i + 1,j,k) + V(i - 1,j,k) + ...
                      + V(i,j - 1,k) + 2*V_Right + V(i,j, k - 1) + ...
                      + 2*V_Back)/hquadro;
    j = 1; k = Nz - 1;
        LapV(i,j,k) = (1/Re)*(-8*V(i,j,k) + V(i + 1,j,k) + V(i - 1,j,k) + ...
                      + V(i,j + 1,k) + 2*V_Left + V(i,j, k - 1) + ...
                      + 2*V_Back)/hquadro;
end

% LapW: 
for k = 2 : Nz - 1
    for j = 2 : Nx - 2
        i = 1;
        LapW(i,j,k) = (1/Re)*(-7*W(i,j,k) + W(i + 1,j,k) + 2*W_Up + ...
                             + W(i,j + 1,k) + W(i,j - 1,k) + W(i,j,k + 1) + ...
                             + W(i,j,k - 1))/hquadro;
        i = Ny - 1; 
        LapW(i,j,k) = (1/Re)*(-7*W(i,j,k) + 2*W_Down + W(i - 1,j,k) + ...
                             + W(i,j + 1,k) + W(i,j - 1,k) + W(i,j,k + 1) + ...
                             + W(i,j,k - 1))/hquadro;
    end
end
for k = 2 : Nz - 1
    for i = 2 : Ny - 2
        j = 1; 
        LapW(i,j,k) = (1/Re)*(-7*W(i,j,k) + W(i + 1,j,k) + W(i - 1,j,k) + ...
                             + W(i,j + 1,k) + 2*W_Left + W(i,j,k + 1) + ...
                             + W(i,j,k - 1))/hquadro;
        j = Nx - 1;
        LapW(i,j,k) = (1/Re)*(-7*W(i,j,k) + W(i + 1,j,k) + W(i + 1,j,k) + ...
                             + 2*W_Right + W(i,j - 1,k) + W(i,j,k + 1) + ...
                             + W(i,j,k - 1))/hquadro;
    end
end

for k = 2 : Nz - 1
    i = 1; j = 1; 
    LapW(i,j,k) = (1/Re)*(-8*W(i,j,k) + W(i + 1,j,k) + 2*W_Left + ...
                         + 2*W_Up + W(i,j + 1,k) + W(i,j,k + 1) + ...
                         + W(i,j,k - 1))/hquadro;
    i = 1; j = Nx - 1;
    LapW(i,j,k) = (1/Re)*(-8*W(i,j,k) + W(i + 1,j,k) + 2*W_Right + ...
                         + W(i,j - 1,k) + 2*W_Up + W(i,j,k + 1) + ...
                         + W(i,j,k - 1))/hquadro;
    i = Ny - 1; j = 1; 
    LapW(i,j,k) = (1/Re)*(-8*W(i,j,k) + 2*W_Left + W(i - 1,j,k) + ...
                         + 2*W_Down + W(i,j + 1,k) + W(i,j,k + 1) + ...
                         + W(i,j,k - 1))/hquadro;
    i = Ny - 1; j = Nx - 1; 
    LapW(i,j,k) = (1/Re)*(-8*W(i,j,k) + 2*W_Right + W(i - 1,j,k) + ...
                         + W(i,j - 1,k) + 2*W_Down + W(i,j,k + 1) + ...
                         + W(i,j,k - 1))/hquadro;     
end

FU = FU + LapU; 
FV = FV + LapV;
FW = FW + LapW;

end

