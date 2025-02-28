clc; clear; close all; 
[Nx,Ny,Nz,x,y,z,h,hquadro,Lx,Ly,Lz,omega] = Params();
[U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
 V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
 W_Up,W_Down,W_Front,W_Back,W_Right,W_Left,Uref]= BCs();
load All_data_Re_100.mat;
load Data_Paper_100.mat;
load Data_Paper_altro_100.mat;

clear title; 

%% Interpolo sui vertici lungo centerline. 
Iu_yz_cen = zeros(Ny-1,Nx-1,Nz-1);
i = 1 : Ny - 2; j = Nx/2; k = 1 : Nz - 1;
Iu_yz_cen(i+1,j,k) = (U(i+1,j,k) + U(i,j,k))/2;

    for k = 1 : Nz - 1
        j = Nx/2; 
        i = Ny;
        Iu_yz_cen(i,j,k) = U_Down;
        i = 1;
        Iu_yz_cen(i,j,k) = U_Up;
    end


Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny; j = Nx/2; k = 1 : Nz - 2;
Iu_ver(i,j,k+1) = (Iu_yz_cen(i,j,k+1) + Iu_yz_cen(i,j,k))/2;

for i = 1 : Ny
    j = Nx/2; 
    k = 1; 
    Iu_ver(i,j,k) = W_Front;
    k = Nz;
    Iu_ver(i,j,k) = W_Back;
end

Iu_ver = squeeze(Iu_ver(:,Nx/2,:));

Iv_xz = zeros(Ny-1,Nx-1,Nz-1); 
i = Ny/2; j = 1 : Nx - 2; k = 1 : Nz - 1; 
Iv_xz(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for k = 1 : Nz - 1
    i = Ny/2;
    j = 1; 
    Iv_xz(i,j,k) = V_Left;
    j = Ny; 
    Iv_xz(i,j,k) = V_Right; 
end

Iv_ver = zeros(Ny,Nx,Nz);
i = Ny/2; j = 1 : Nx; k = 1 : Nz - 2; 
Iv_ver(i,j,k + 1) = (Iv_xz(i,j,k + 1) + Iv_xz(i,j,k))/2;

for j = 1 : Nx 
    i = Ny/2; 
    k = 1; 
    Iv_ver(i,j,k) = V_Front; 
    k = Nz; 
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(Ny/2,:,:));

%% Plotto centerline. 
x_paper = Data_Paper_100(:,1); 
u_paper = Data_Paper_100(:,2);

y_paper = Data_paper_altro_100(:,1);
v_paper = Data_paper_altro_100(:,2);

centerline_u = squeeze(Iu_ver(:,Nz/2));
centerline_v = squeeze(Iv_ver(:,Nz/2));

% Figura 1: U at (y, 0.5, 0.5)
figure(1);
plot(x(1:Ny), centerline_u, '-k', 'LineWidth', 1.2); hold on;  % Linea nera con spessore moderato
plot(x_paper, u_paper * 2 - 1, 'o', 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 6);  % Punti grigi chiari

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità U', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di U al centro del dominio (y, 0.5, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva U', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara

saveas(gcf, 'Figura_U_centerline.png');

% Figura 2: V at (0.5, x, 0.5)
figure(2);
plot(x(1:Ny), centerline_v, '-', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.2); hold on;  % Linea blu tenue
plot(y_paper + 0.5, v_paper * 2, 's', 'MarkerEdgeColor', [0.3 0.5 0.4], 'MarkerFaceColor', [0.6 0.8 0.7], 'MarkerSize', 6);  % Punti verde tenue

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità V', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di V al centro del dominio (0.5, x, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva V', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara

saveas(gcf, 'Figura_V_centerline.png');

figcont = 3;
%% Interpolo sui vertici ad x = per il quiver.
figcont = 3;
valori = [0.962,0.854,0.5,0.1446];
for cont = 1 : 4
val = valori(cont);

Iw_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = round(val*(Nz-1)); k = 1 : Nz;
Iw_ver(i+1,j,k) = (W(i + 1,j,k) + W(i,j,k))/2;

for k = 1 : Nz
    j = round(val*(Nz-1)); 
    i = 1; 
    Iw_ver(i,j,k) = W_Up;
    i = Ny;
    Iw_ver(i,j,k) = W_Down;
end

Iw_ver = squeeze(Iw_ver(:,round(val*(Nz-1)),:));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = round(val*(Nz-1)); k = 1 : Nz - 2; 
Iv_ver(i,j,k+1) = (V(i,j,k+1) + V(i,j,k))/2;

for i = 1 : Ny 
    j = round(val*(Nz-1));
    k = 1; 
    Iv_ver(i,j,k) = V_Front;
    k = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,round(val*(Nz-1)),:));

% Crea una nuova figura
figure(figcont); 
figcont = figcont + 1;

% Meshgrid per i dati
[Z, Y] = meshgrid(z, y);

% Titolo del grafico
titolo = strcat('Quiver at x = ', num2str(val));

% Grafico del campo vettoriale
quiver(Z, Y, Iw_ver, Iv_ver, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'MaxHeadSize', 0.3);  % Frecce grigie morbide

% Etichette degli assi e titolo
xlabel('Coordinata Z', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coordinata Y', 'FontSize', 12, 'FontWeight', 'bold');
title(titolo, 'FontSize', 14, 'FontWeight', 'bold');

% Miglioramenti estetici
grid on;  % Attiva la griglia
set(gca, 'FontSize', 10);  % Dimensione font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara
axis tight;  % Adatta gli assi ai dati

nome_file = strcat('Quiver_x_', num2str(val), '.png');
saveas(gcf, nome_file);  % Salva la figura come immagine PNG

end


%% Interpolo a z = per il plot 3D. 

figure(figcont);
set(gcf, 'Color', [0.3, 0.3, 0.3]);

figcont = figcont + 1;
valori = [0.1,0.5,0.90];
for cont = 1 : 3
val = valori(cont);

Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = 1 : Nx; k = round(val*(Nz-1));
Iu_ver(i + 1,j,k) = (U(i + 1,j,k) + U(i,j,k))/2;

for j = 1 : Nx
    k = round(val*(Nz-1)); 
    i = 1; 
    Iu_ver(i,j,k) = U_Up;
    i = Ny;
    Iu_ver(i,j,k) = U_Down;
end

Iu_ver = squeeze(Iu_ver(:,:,round(val*(Nz-1))));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = 1 : Nx - 2; k = round(val*(Nz-1)); 
Iv_ver(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for i = 1 : Ny 
    k = round(val*(Nz-1));
    j = 1; 
    Iv_ver(i,j,k) = V_Front;
    j = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,:,round(val*(Nz-1))));

% Creazione della figura e delle streamline

[X, Y] = meshgrid(x, y);
a = streamslice(X, Y, Iu_ver, Iv_ver, 'arrows', 'none');
hold on;

% Colormap per colorare le streamline
cmap = turbo(128);  % Colormap simile al gradiente fluido della figura

% Calcola la velocità nel campo
velocity_magnitude = sqrt(Iu_ver.^2 + Iv_ver.^2);

% Trova i limiti delle velocità per la normalizzazione
min_val = min(velocity_magnitude(:));
max_val = max(velocity_magnitude(:));

% Itera sulle streamline
for k = 1:length(a)
    % Estrai le coordinate della streamline originale
    xData = a(k).XData;
    yData = a(k).YData;
    zData = val * ones(size(xData));  % Streamline inizialmente nel piano X-Y
    
    % Interpola il valore della velocità per ogni punto della streamline
    velocity_values = interp2(X, Y, velocity_magnitude, xData, yData, 'linear', NaN);
    
    % Normalizza i valori di velocità per mappare nella colormap
    normalized_velocity = (velocity_values - min_val) / (max_val - min_val);
    normalized_velocity(isnan(normalized_velocity)) = 0;  % Gestione NaN
    color_indices = round(normalized_velocity * (size(cmap, 1) - 1)) + 1;
    
    % Assicurati che gli indici siano validi
    color_indices = min(max(color_indices, 1), size(cmap, 1));
    
    % Disegna la streamline con il gradiente
    for j = 1:length(xData) - 1
        % Colore corrente dalla colormap
        current_color = cmap(color_indices(j), :);
        
        % Disegna il segmento della streamline
        plot3([xData(j), xData(j+1)], ...
              [yData(j), yData(j+1)], ...
              [zData(j), zData(j+1)], ...
              'Color', current_color, 'LineWidth', 1.5);
    end
end

% Configurazione degli assi e sfondo colorato
ax = gca;
ax.Color = [0.3, 0.3, 0.3];  % Colore di sfondo degli assi
ax.XColor = 'w';  % Colore dei numeri sugli assi
ax.YColor = 'w';
ax.ZColor = 'w';
xlabel('X', 'FontWeight', 'bold', 'Color', 'w');
ylabel('Y', 'FontWeight', 'bold', 'Color', 'w');
zlabel('Z', 'FontWeight', 'bold', 'Color', 'w');
zlim([0.01, 1]);
xlim([0, 1]);
ylim([0, 1]);
grid on;
set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.3);  % Griglia chiara ma semitrasparente

% Disegna il cubo attorno alle streamline
x_limits = xlim;
y_limits = ylim;
z_limits = zlim;

% Linee orizzontali (bordi inferiori e superiori del cubo)
for z = [z_limits(1), z_limits(2)]
    line([x_limits(1), x_limits(2)], [y_limits(1), y_limits(1)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(2)], [y_limits(2), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(1)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(2), x_limits(2)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
end

% Linee verticali
for xx = [x_limits(1), x_limits(2)]
    for yy = [y_limits(1), y_limits(2)]
        line([xx, xx], [yy, yy], [z_limits(1), z_limits(2)], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    end
end

% Barra di colore per indicare la scala delle velocità
c = colorbar;
c.Label.String = 'Velocity';
c.Label.FontWeight = 'bold';
c.Label.Color = 'w';
set(c, 'Color', 'w');  % Colore bianco per la barra di colore

% Vista 3D elegante
view(3);  
hold off;



end
%% Cambio il Re; 

clearvars -except figcont; close all;
[Nx,Ny,Nz,x,y,z,h,hquadro,Lx,Ly,Lz,omega] = Params();
[U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
 V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
 W_Up,W_Down,W_Front,W_Back,W_Right,W_Left,Uref]= BCs();
load All_data_Re_400.mat;
load Data_Paper_400.mat;
load Data_Paper_altro_400.mat;

clear title;

%% Interpolo sui vertici. 
Iu_yz_cen = zeros(Ny-1,Nx-1,Nz-1);
i = 1 : Ny - 2; j = Nx/2; k = 1 : Nz - 1;
Iu_yz_cen(i+1,j,k) = (U(i+1,j,k) + U(i,j,k))/2;

    for k = 1 : Nz - 1
        j = Nx/2; 
        i = Ny;
        Iu_yz_cen(i,j,k) = U_Down;
        i = 1;
        Iu_yz_cen(i,j,k) = U_Up;
    end


Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny; j = Nx/2; k = 1 : Nz - 2;
Iu_ver(i,j,k+1) = (Iu_yz_cen(i,j,k+1) + Iu_yz_cen(i,j,k))/2;

for i = 1 : Ny
    j = Nx/2; 
    k = 1; 
    Iu_ver(i,j,k) = U_Front;
    k = Nz;
    Iu_ver(i,j,k) = U_Back;
end

Iu_ver = squeeze(Iu_ver(:,Nx/2,:));

Iv_xz = zeros(Ny-1,Nx-1,Nz-1); 
i = Ny/2; j = 1 : Nx - 2; k = 1 : Nz - 1; 
Iv_xz(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for k = 1 : Nz - 1
    i = Ny/2;
    j = 1; 
    Iv_xz(i,j,k) = V_Left;
    j = Ny; 
    Iv_xz(i,j,k) = V_Right; 
end

Iv_ver = zeros(Ny,Nx,Nz);
i = Ny/2; j = 1 : Nx; k = 1 : Nz - 2; 
Iv_ver(i,j,k + 1) = (Iv_xz(i,j,k + 1) + Iv_xz(i,j,k))/2;

for j = 1 : Nx 
    i = Ny/2; 
    k = 1; 
    Iv_ver(i,j,k) = V_Front; 
    k = Nz; 
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(Ny/2,:,:));

%% Plotto le centerline. 
x_paper = Data_Paper_400(:,1); 
u_paper = Data_Paper_400(:,2);

y_paper = Data_paper_altro_400(:,1);
v_paper = Data_paper_altro_400(:,2);

centerline_u = squeeze(Iu_ver(:,Nz/2));
centerline_v = squeeze(Iv_ver(:,Nz/2));

% Figura 1: U at (y, 0.5, 0.5)
figure(figcont);
figcont = figcont + 1;
plot(x(1:Ny), centerline_u, '-k', 'LineWidth', 1.2); hold on;  % Linea nera con spessore moderato
plot(x_paper, u_paper * 2 - 1, 'o', 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 6);  % Punti grigi chiari

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità U', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di U al centro del dominio (y, 0.5, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva U', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara

saveas(gcf, 'Figura_U_centerline.png');

% Figura 2: V at (0.5, x, 0.5)
figure(figcont);
figcont = figcont + 1;
plot(x(1:Ny), centerline_v, '-', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.2); hold on;  % Linea blu tenue
plot(y_paper + 0.5, v_paper * 2, 's', 'MarkerEdgeColor', [0.3 0.5 0.4], 'MarkerFaceColor', [0.6 0.8 0.7], 'MarkerSize', 6);  % Punti verde tenue

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità V', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di V al centro del dominio (0.5, x, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva V', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara

saveas(gcf, 'Figura_V_centerline.png');


%% Interpolo sui vertici ad x = per il quiver.
valori = [0.962,0.854,0.5,0.1446];
for cont = 1 : 4
val = valori(cont);

Iw_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = round(val*(Nz-1)); k = 1 : Nz;
Iw_ver(i+1,j,k) = (W(i + 1,j,k) + W(i,j,k))/2;

for k = 1 : Nz
    j = round(val*(Nz-1)); 
    i = 1; 
    Iw_ver(i,j,k) = W_Up;
    i = Ny;
    Iw_ver(i,j,k) = W_Down;
end

Iw_ver = squeeze(Iw_ver(:,round(val*(Nz-1)),:));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = round(val*(Nz-1)); k = 1 : Nz - 2; 
Iv_ver(i,j,k+1) = (V(i,j,k+1) + V(i,j,k))/2;

for i = 1 : Ny 
    j = round(val*(Nz-1));
    k = 1; 
    Iv_ver(i,j,k) = V_Front;
    k = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,round(val*(Nz-1)),:));

% Crea una nuova figura
figure(figcont); 
figcont = figcont + 1;

% Meshgrid per i dati
[Z, Y] = meshgrid(z, y);

% Titolo del grafico
titolo = strcat('Quiver at x = ', num2str(val));

% Grafico del campo vettoriale
quiver(Z, Y, Iw_ver, Iv_ver, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'MaxHeadSize', 0.3);  % Frecce grigie morbide

% Etichette degli assi e titolo
xlabel('Coordinata Z', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coordinata Y', 'FontSize', 12, 'FontWeight', 'bold');
title(titolo, 'FontSize', 14, 'FontWeight', 'bold');

% Miglioramenti estetici
grid on;  % Attiva la griglia
set(gca, 'FontSize', 10);  % Dimensione font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara
axis tight;  % Adatta gli assi ai dati

nome_file = strcat('Quiver_x_', num2str(val), '.png');
saveas(gcf, nome_file);  % Salva la figura come immagine PNG

end

%% Interpolo a z = per il plot 3D. 
figure(figcont); 
set(gcf, 'Color', [0.3, 0.3, 0.3]);
figcont = figcont + 1;
valori = [0.1,0.5,0.90];
for cont = 1 : 3
val = valori(cont);

Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = 1 : Nx; k = round(val*(Nz-1));
Iu_ver(i + 1,j,k) = (U(i + 1,j,k) + U(i,j,k))/2;

for j = 1 : Nx
    k = round(val*(Nz-1)); 
    i = 1; 
    Iu_ver(i,j,k) = U_Up;
    i = Ny;
    Iu_ver(i,j,k) = U_Down;
end

Iu_ver = squeeze(Iu_ver(:,:,round(val*(Nz-1))));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = 1 : Nx - 2; k = round(val*(Nz-1)); 
Iv_ver(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for i = 1 : Ny 
    k = round(val*(Nz-1));
    j = 1; 
    Iv_ver(i,j,k) = V_Front;
    j = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,:,round(val*(Nz-1))));

% Creazione della figura e delle streamline

[X, Y] = meshgrid(x, y);
a = streamslice(X, Y, Iu_ver, Iv_ver, 'arrows', 'none');
hold on;

% Colormap per colorare le streamline
cmap = turbo(128);  % Colormap simile al gradiente fluido della figura

% Calcola la velocità nel campo
velocity_magnitude = sqrt(Iu_ver.^2 + Iv_ver.^2);

% Trova i limiti delle velocità per la normalizzazione
min_val = min(velocity_magnitude(:));
max_val = max(velocity_magnitude(:));

% Itera sulle streamline
for k = 1:length(a)
    % Estrai le coordinate della streamline originale
    xData = a(k).XData;
    yData = a(k).YData;
    zData = val * ones(size(xData));  % Streamline inizialmente nel piano X-Y
    
    % Interpola il valore della velocità per ogni punto della streamline
    velocity_values = interp2(X, Y, velocity_magnitude, xData, yData, 'linear', NaN);
    
    % Normalizza i valori di velocità per mappare nella colormap
    normalized_velocity = (velocity_values - min_val) / (max_val - min_val);
    normalized_velocity(isnan(normalized_velocity)) = 0;  % Gestione NaN
    color_indices = round(normalized_velocity * (size(cmap, 1) - 1)) + 1;
    
    % Assicurati che gli indici siano validi
    color_indices = min(max(color_indices, 1), size(cmap, 1));
    
    % Disegna la streamline con il gradiente
    for j = 1:length(xData) - 1
        % Colore corrente dalla colormap
        current_color = cmap(color_indices(j), :);
        
        % Disegna il segmento della streamline
        plot3([xData(j), xData(j+1)], ...
              [yData(j), yData(j+1)], ...
              [zData(j), zData(j+1)], ...
              'Color', current_color, 'LineWidth', 1.5);
    end
end

% Configurazione degli assi e sfondo colorato
ax = gca;
ax.Color = [0.3, 0.3, 0.3];  % Colore di sfondo degli assi
ax.XColor = 'w';  % Colore dei numeri sugli assi
ax.YColor = 'w';
ax.ZColor = 'w';
xlabel('X', 'FontWeight', 'bold', 'Color', 'w');
ylabel('Y', 'FontWeight', 'bold', 'Color', 'w');
zlabel('Z', 'FontWeight', 'bold', 'Color', 'w');
zlim([0.01, 1]);
xlim([0, 1]);
ylim([0, 1]);
grid on;
set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.3);  % Griglia chiara ma semitrasparente

% Disegna il cubo attorno alle streamline
x_limits = xlim;
y_limits = ylim;
z_limits = zlim;

% Linee orizzontali (bordi inferiori e superiori del cubo)
for z = [z_limits(1), z_limits(2)]
    line([x_limits(1), x_limits(2)], [y_limits(1), y_limits(1)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(2)], [y_limits(2), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(1)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(2), x_limits(2)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
end

% Linee verticali
for xx = [x_limits(1), x_limits(2)]
    for yy = [y_limits(1), y_limits(2)]
        line([xx, xx], [yy, yy], [z_limits(1), z_limits(2)], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    end
end

% Barra di colore per indicare la scala delle velocità
c = colorbar;
c.Label.String = 'Velocity';
c.Label.FontWeight = 'bold';
c.Label.Color = 'w';
set(c, 'Color', 'w');  % Colore bianco per la barra di colore

% Vista 3D elegante
view(3);  
hold off;




end

%% Cambio il Re; 
clearvars -except figcont; close all;  
[Nx,Ny,Nz,x,y,z,h,hquadro,Lx,Ly,Lz,omega] = Params();
[U_Up,U_Down,U_Front,U_Back,U_Right,U_Left, ...
 V_Up,V_Down,V_Front,V_Back,V_Right,V_Left, ...
 W_Up,W_Down,W_Front,W_Back,W_Right,W_Left,Uref]= BCs();
load All_data_Re_1000.mat;
load Data_Paper_1000.mat;
load Data_Paper_altro_1000.mat;

clear title; 
%% Interpolo sui vertici. 
Iu_yz_cen = zeros(Ny-1,Nx-1,Nz-1);
i = 1 : Ny - 2; j = Nx/2; k = 1 : Nz - 1;
Iu_yz_cen(i+1,j,k) = (U(i+1,j,k) + U(i,j,k))/2;

    for k = 1 : Nz - 1
        j = Nx/2; 
        i = Ny;
        Iu_yz_cen(i,j,k) = U_Down;
        i = 1;
        Iu_yz_cen(i,j,k) = U_Up;
    end


Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny; j = Nx/2; k = 1 : Nz - 2;
Iu_ver(i,j,k+1) = (Iu_yz_cen(i,j,k+1) + Iu_yz_cen(i,j,k))/2;

for i = 1 : Ny
    j = Nx/2; 
    k = 1; 
    Iu_ver(i,j,k) = U_Front;
    k = Nz;
    Iu_ver(i,j,k) = U_Back;
end

Iu_ver = squeeze(Iu_ver(:,Nx/2,:));

Iv_xz = zeros(Ny-1,Nx-1,Nz-1); 
i = Ny/2; j = 1 : Nx - 2; k = 1 : Nz - 1; 
Iv_xz(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for k = 1 : Nz - 1
    i = Ny/2;
    j = 1; 
    Iv_xz(i,j,k) = V_Left;
    j = Ny; 
    Iv_xz(i,j,k) = V_Right; 
end

Iv_ver = zeros(Ny,Nx,Nz);
i = Ny/2; j = 1 : Nx; k = 1 : Nz - 2; 
Iv_ver(i,j,k + 1) = (Iv_xz(i,j,k + 1) + Iv_xz(i,j,k))/2;

for j = 1 : Nx 
    i = Ny/2; 
    k = 1; 
    Iv_ver(i,j,k) = V_Front; 
    k = Nz; 
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(Ny/2,:,:));

%% Plotto le centerline. 
x_paper = Data_paper_Re_1000(:,1); 
u_paper = Data_paper_Re_1000(:,2);

y_paper = Data_paper_altro_1000(:,1);
v_paper = Data_paper_altro_1000(:,2);

centerline_u = squeeze(Iu_ver(:,Nz/2));
centerline_v = squeeze(Iv_ver(:,Nz/2));

% Figura 1: U at (y, 0.5, 0.5)
figure(figcont);
figcont = figcont + 1;
plot(x(1:Ny), centerline_u, '-k', 'LineWidth', 1.2); hold on;  % Linea nera con spessore moderato
plot(x_paper, u_paper * 2 - 1, 'o', 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 6);  % Punti grigi chiari

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità U', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di U al centro del dominio (y, 0.5, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva U', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara


saveas(gcf, 'Figura_U_centerline.png');

% Figura 2: V at (0.5, x, 0.5)
figure(figcont);
figcont = figcont + 1;
axis([0 1 0 1]);
plot(x(1:Ny), centerline_v, '-', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.2); hold on;  % Linea blu tenue
plot(y_paper + 0.5, v_paper * 2, 's', 'MarkerEdgeColor', [0.3 0.5 0.4], 'MarkerFaceColor', [0.6 0.8 0.7], 'MarkerSize', 6);  % Punti verde tenue

grid on;  % Attiva la griglia
xlabel('Asse X', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Velocità V', 'FontSize', 12, 'FontWeight', 'bold'); 
title('Distribuzione di V al centro del dominio (0.5, x, 0.5)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Curva V', 'Dati sperimentali', 'Location', 'best');  % Legenda automatica
set(gca, 'FontSize', 10);  % Aumenta il font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara



saveas(gcf, 'Figura_V_centerline.png');




%% Interpolo sui vertici ad x = per il quiver.
valori = [0.962,0.854,0.5,0.1446];
for cont = 1 : 4
val = valori(cont);

Iw_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = round(val*(Nz-1)); k = 1 : Nz;
Iw_ver(i+1,j,k) = (W(i + 1,j,k) + W(i,j,k))/2;

for k = 1 : Nz
    j = round(val*(Nz-1)); 
    i = 1; 
    Iw_ver(i,j,k) = W_Up;
    i = Ny;
    Iw_ver(i,j,k) = W_Down;
end

Iw_ver = squeeze(Iw_ver(:,round(val*(Nz-1)),:));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = round(val*(Nz-1)); k = 1 : Nz - 2; 
Iv_ver(i,j,k+1) = (V(i,j,k+1) + V(i,j,k))/2;

for i = 1 : Ny 
    j = round(val*(Nz-1));
    k = 1; 
    Iv_ver(i,j,k) = V_Front;
    k = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,round(val*(Nz-1)),:));

% Crea una nuova figura
figure(figcont); 
figcont = figcont + 1;

% Meshgrid per i dati
[Z, Y] = meshgrid(z, y);

% Titolo del grafico
titolo = strcat('Quiver at x = ', num2str(val));

% Grafico del campo vettoriale
quiver(Z, Y, Iw_ver, Iv_ver, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'MaxHeadSize', 0.3);  % Frecce grigie morbide

% Etichette degli assi e titolo
xlabel('Coordinata Z', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coordinata Y', 'FontSize', 12, 'FontWeight', 'bold');
title(titolo, 'FontSize', 14, 'FontWeight', 'bold');

% Miglioramenti estetici
grid on;  % Attiva la griglia
set(gca, 'FontSize', 10);  % Dimensione font degli assi
set(gca, 'GridColor', [0.8 0.8 0.8]);  % Griglia grigia chiara
axis tight;  % Adatta gli assi ai dati

nome_file = strcat('Quiver_x_', num2str(val), '.png');
saveas(gcf, nome_file);  % Salva la figura come immagine PNG

end

%% Interpolo a z = per il plot 3D. 
figure(figcont); 
set(gcf, 'Color', [0.3, 0.3, 0.3]);
figcont = figcont + 1;
valori = [0.1,0.5,0.90];
for cont = 1 : 3
val = valori(cont);

Iu_ver = zeros(Ny,Nx,Nz);
i = 1 : Ny-2; j = 1 : Nx; k = round(val*(Nz-1));
Iu_ver(i + 1,j,k) = (U(i + 1,j,k) + U(i,j,k))/2;

for j = 1 : Nx
    k = round(val*(Nz-1)); 
    i = 1; 
    Iu_ver(i,j,k) = U_Up;
    i = Ny;
    Iu_ver(i,j,k) = U_Down;
end

Iu_ver = squeeze(Iu_ver(:,:,round(val*(Nz-1))));

Iv_ver = zeros(Ny,Nx-1,Nz-1);
i = 1 : Ny; j = 1 : Nx - 2; k = round(val*(Nz-1)); 
Iv_ver(i,j+1,k) = (V(i,j+1,k) + V(i,j,k))/2;

for i = 1 : Ny 
    k = round(val*(Nz-1));
    j = 1; 
    Iv_ver(i,j,k) = V_Front;
    j = Nz;
    Iv_ver(i,j,k) = V_Back;
end

Iv_ver = squeeze(Iv_ver(:,:,round(val*(Nz-1))));

% Creazione della figura e delle streamline

[X, Y] = meshgrid(x, y);
a = streamslice(X, Y, Iu_ver, Iv_ver, 'arrows', 'none');
hold on;

% Colormap per colorare le streamline
cmap = turbo(128);  % Colormap simile al gradiente fluido della figura

% Calcola la velocità nel campo
velocity_magnitude = sqrt(Iu_ver.^2 + Iv_ver.^2);

% Trova i limiti delle velocità per la normalizzazione
min_val = min(velocity_magnitude(:));
max_val = max(velocity_magnitude(:));

% Itera sulle streamline
for k = 1:length(a)
    % Estrai le coordinate della streamline originale
    xData = a(k).XData;
    yData = a(k).YData;
    zData = val * ones(size(xData));  % Streamline inizialmente nel piano X-Y
    
    % Interpola il valore della velocità per ogni punto della streamline
    velocity_values = interp2(X, Y, velocity_magnitude, xData, yData, 'linear', NaN);
    
    % Normalizza i valori di velocità per mappare nella colormap
    normalized_velocity = (velocity_values - min_val) / (max_val - min_val);
    normalized_velocity(isnan(normalized_velocity)) = 0;  % Gestione NaN
    color_indices = round(normalized_velocity * (size(cmap, 1) - 1)) + 1;
    
    % Assicurati che gli indici siano validi
    color_indices = min(max(color_indices, 1), size(cmap, 1));
    
    % Disegna la streamline con il gradiente
    for j = 1:length(xData) - 1
        % Colore corrente dalla colormap
        current_color = cmap(color_indices(j), :);
        
        % Disegna il segmento della streamline
        plot3([xData(j), xData(j+1)], ...
              [yData(j), yData(j+1)], ...
              [zData(j), zData(j+1)], ...
              'Color', current_color, 'LineWidth', 1.5);
    end
end

% Configurazione degli assi e sfondo colorato
ax = gca;
ax.Color = [0.3, 0.3, 0.3];  % Colore di sfondo degli assi
ax.XColor = 'w';  % Colore dei numeri sugli assi
ax.YColor = 'w';
ax.ZColor = 'w';
xlabel('X', 'FontWeight', 'bold', 'Color', 'w');
ylabel('Y', 'FontWeight', 'bold', 'Color', 'w');
zlabel('Z', 'FontWeight', 'bold', 'Color', 'w');
zlim([0.01, 1]);
xlim([0, 1]);
ylim([0, 1]);
grid on;
set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.3);  % Griglia chiara ma semitrasparente

% Disegna il cubo attorno alle streamline
x_limits = xlim;
y_limits = ylim;
z_limits = zlim;

% Linee orizzontali (bordi inferiori e superiori del cubo)
for z = [z_limits(1), z_limits(2)]
    line([x_limits(1), x_limits(2)], [y_limits(1), y_limits(1)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(2)], [y_limits(2), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(1), x_limits(1)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    line([x_limits(2), x_limits(2)], [y_limits(1), y_limits(2)], [z, z], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
end

% Linee verticali
for xx = [x_limits(1), x_limits(2)]
    for yy = [y_limits(1), y_limits(2)]
        line([xx, xx], [yy, yy], [z_limits(1), z_limits(2)], 'Color', 'w', 'LineStyle', '-', 'LineWidth', 1.5);
    end
end

% Barra di colore per indicare la scala delle velocità
c = colorbar;
c.Label.String = 'Velocity';
c.Label.FontWeight = 'bold';
c.Label.Color = 'w';
set(c, 'Color', 'w');  % Colore bianco per la barra di colore

% Vista 3D elegante
view(3);  
hold off;


end