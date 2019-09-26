function ResultPlots = ComsolCompareHoles(Nk, Nmesh, lbda, Mpoint, Kpoint)
%% ComsolCompare function creates plots of COMSOL results and FDFD results for comparison
%   Different k space ranges are plotted on different plots in order to
%   correspond to the different incident angles as simulated in COMSOL.

%% Inputs
%Nk = Number of total k points taken from the band structure main code
%Nmesh = Number of total mesh points from the band structure main code
%lbda = (Nk*Nmesh) Matrix of wavelengths from the eigenstates calculated in main code (nm)
%MPoint = k-coordinates of the high symmetry point in the irreducible BZ
%KPoint = k-coordinates of the high symmetry point in the irreducible BZ
FEMdata0 = importdata('FEMdataHoles0.mat'); % Matrix of wavelength and reflectance COMSOL values for 0º incidence
FEMdata30 = importdata('FEMdataHoles30.mat'); % Matrix of wavelength and reflectance COMSOL values for 45º incidence

%% Gamma to M (k-space) or 0º Incidence (real space)
% Plot band structure as k vs wavelength
subplot(2,2,1); % Sets up 2 rows, 2 columns of plots in a single figure
kPlotVector = zeros(Nk, 1); % Initializes y-axis (k)
for i = 1:Nk % Cycle through all k points
    kPlotVector(i) = i; % Counting k points
end

PlotModes = zeros(Nk, 1); % Initialize plot?
for i = 1:Nmesh % Cycle through number of mesh points
    PlotModes(i) = plot(lbda(:, i), kPlotVector, 'b'); %k vs wavelength
    hold on;
end

set(gca, 'Ytick', [1, Mpoint]); % Create y-axis labels
set(gca, 'YTickLabel', ['G'; 'M']);
axis([min(FEMdata0(:, 1)) max(FEMdata0(:, 1)) 1 Mpoint]); % Set y-axis range
title('Dispersion Relation Gamma to M');
ylabel('k');
xlabel('Wavelength (nm)');

% Plot FEM data as reflectance vs wavelength
subplot(2,2,3);
plot(FEMdata0(:, 1), FEMdata0(:, 2));
axis([min(FEMdata0(:, 1)) max(FEMdata0(:, 1)) 0 1]);
title('Reflectance vs Wavelength at 0º');
ylabel('Reflectance');
xlabel('Wavelength (nm)');

%% K to Gamma (k-space) or 30º Incidence (real space)
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

set(gca, 'Ytick', [Kpoint, Nk]); % Create y-axis labels
set(gca, 'YTickLabel', ['K'; 'G']);
axis([min(FEMdata30(:, 1)) max(FEMdata30(:, 1)) Kpoint Nk]) % Set y-axis range
title('Dispersion Relation K to Gamma');
ylabel('k');
xlabel('Wavelength (nm)');

% Plot FEM data as reflectance vs wavelength
subplot(2,2,4);
plot(FEMdata30(:, 1), FEMdata30(:, 2));
axis([min(FEMdata30(:, 1)) max(FEMdata30(:, 1)) 0 1])
title('Reflectance vs Wavelength at 30º');
ylabel('Reflectance');
xlabel('Wavelength (nm)');

end