function ResultPlots = ComsolComparePosts(Nk, Nmesh, lbda, Xpoint)
%% ComsolCompare function creates plots of COMSOL results and FDFD results for comparison
%   Different k space ranges are plotted on different plots in order to
%   correspond to the different incident angles as simulated in COMSOL.

%% Inputs
%Nk = Number of total k points taken from the band structure main code
%Nmesh = Number of total mesh points from the band structure main code
%lbda = (Nk*Nmesh) Matrix of wavelengths from the eigenstates calculated in main code (nm)
%XPoint = k-coordinates of the first high symmetry point in the irreducible BZ
FEMdata0 = importdata('FEMdataPosts0.mat'); % Matrix of wavelength and reflectance COMSOL values for 0º incidence
FEMdata45 = importdata('FEMdataPosts45.mat'); % Matrix of wavelength and reflectance COMSOL values for 45º incidence

%% Gamma to X (k-space) or 0º Incidence (real space)
% Plot band structure as k vs wavelength
subplot(2,2,1); % Sets up 2 rows, 1 column of plots in a single figure
kPlotVector = zeros(Nk, 1); % Initializes y-axis (k)
for i = 1:Nk % Cycle through all k points
    kPlotVector(i) = i; % Counting k points
end

PlotModes = zeros(Nk, 1); % Initialize plot?
for i = 1:Nmesh % Cycle through number of mesh points
    PlotModes(i) = plot(lbda(:, i), kPlotVector, 'b'); %k vs wavelength
    hold on;
end

set(gca, 'Ytick', [1, Xpoint]); % Create y-axis labels
set(gca, 'YTickLabel', ['G'; 'X']);
axis([min(FEMdata0(:, 1)) max(FEMdata0(:, 1)) 1 Xpoint]); % Set y-axis range
title('Dispersion Relation Gamma to X');
ylabel('k');
xlabel('Wavelength (nm)');

% Plot FEM data as reflectance vs wavelength
subplot(2,2,3);
plot(FEMdata0(:, 1), FEMdata0(:, 2));
axis([min(FEMdata0(:, 1)) max(FEMdata0(:, 1)) 0 1]);
title('Reflectance vs Wavelength at 0º');
ylabel('Reflectance');
xlabel('Wavelength (nm)');

%% S/M to Gamma (k-space) or 45º Incidence (real space)
% Plot band structure as k vs wavelength
subplot(2,2,2); % Sets up 2 rows, 1 column of plots in a single figure
kPlotVector = zeros(Nk, 1); % Initializes y-axis (k)
for i = 1:Nk % Cycle through all k points
    kPlotVector(i) = i; % Counting k points
end

PlotModes = zeros(Nk, 1); % Initialize plot?
for i = 1:Nmesh % Cycle through number of mesh points
    PlotModes(i) = plot(lbda(:, i), kPlotVector, 'b'); %k vs wavelength
    hold on;
end

Spoint = 2 * Xpoint;
set(gca, 'Ytick', [Spoint, Nk]); % Create y-axis labels
set(gca, 'YTickLabel', ['S'; 'G']);
axis([min(FEMdata45(:, 1)) max(FEMdata45(:, 1)) Spoint Nk]) % Set y-axis range
title('Dispersion Relation S/M to Gamma');
ylabel('k');
xlabel('Wavelength (nm)');

% Plot FEM data as reflectance vs wavelength
subplot(2,2,4);
plot(FEMdata45(:, 1), FEMdata45(:, 2));
axis([min(FEMdata45(:, 1)) max(FEMdata45(:, 1)) 0 1])
title('Reflectance vs Wavelength at 45º');
ylabel('Reflectance');
xlabel('Wavelength (nm)');

end