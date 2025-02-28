function [Px,Py,Pz] = Grad_Calc(P)
[Nx,Ny,Nz,~,~,~,h] = Params();

Px = zeros(Ny - 1, Nx, Nz - 1);  
Py = zeros(Ny, Nx - 1, Nz - 1);
Pz = zeros(Ny - 1, Nx - 1, Nz);

for i = 1 : Ny - 1
    for j = 2 : Nx - 1
        for k = 1 : Nz - 1
            Px(i,j,k) = (P(i,j,k) - P(i,j - 1,k))/h; 
        end
    end
end

for i = 2 : Ny - 1
    for j = 1 : Nx - 1
        for k = 1 : Nz - 1
            Py(i,j,k) = (P(i,j,k) - P(i - 1,j,k))/h;
        end
    end
end

for i = 1 : Ny - 1
    for j = 1 : Nx - 1
        for k = 2 : Nz - 1
            Pz(i,j,k) = (P(i,j,k) - P(i,j,k - 1))/h;
        end
    end
end

end

