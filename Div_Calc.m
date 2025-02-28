function [Div] = Div_Calc(U,V,W)
[Nx,Ny,Nz,~,~,~,h] = Params();

Div = zeros(Ny - 1, Nx - 1, Nz - 1); 

for i = 1 : Ny - 1
    for j = 1 : Nx - 1
        for k = 1 : Nz - 1
            Div(i,j,k) = (U(i,j + 1,k) - U(i,j,k))/h + ...
                       + (V(i + 1,j,k) - V(i,j,k))/h + ...
                       + (W(i,j,k + 1) - W(i,j,k))/h;
        end
    end
end

end

