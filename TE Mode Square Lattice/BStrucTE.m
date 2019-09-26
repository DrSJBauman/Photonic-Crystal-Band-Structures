clear all
clc
clf

%% Parameters:  
Period = 300; % Distance between atoms [m]
P1 = Period; % Distance along basis vector 1 [m]
P2 = Period; % Distance along basis vector 2 [m]
DPoints = 25; % Number of k points along the diagonal of the calculated irreducible BZ
sqPix = 15; %  Number of mesh points in the photonic crystal unit cell in one direction
xPix = sqPix; % Number of mesh points along x
yPix = sqPix; % Number of mesh points along y
T = 2;  % How many Pix^2 matrices concatenated in one row?
wmax = 0.8; % Maximum y-axis value on the final band structure plot
EpsC = 8.9; % Relative dielectric constant of the cylinders (posts/holes) [unitless]
EpsB = 1; % Relative dielectric constant of the bulk surrounding the posts/holes [unitless]
theta = 90; % Angle between the real space basis vectors [degrees]
rRatio = 0.2; % Ratio between the radius of the cylinders and the period  
if rRatio > 0.5
    error('Settle down there cowboy. Try a smaller rRatio.');
end
c =  2.99792458 * 10^8; % Speed of light [m/s]
Nmesh = xPix * yPix; % Number of total mesh elements in the unit cell in real space

%% Create unit cell by applying dielectric constant values to mesh positions.
UnitCell = GenGeom(P1, P2, theta, xPix, yPix, EpsC, EpsB, rRatio);
% UnitCell = ones(16, 3);
% UnitCell(6, 3) = EpsC;
% UnitCell(7, 3) = EpsC;
% UnitCell(10, 3) = EpsC;
% UnitCell(11, 3) = EpsC;
% Testshe 
%UnitCell(:,3)=1;
%UnitCell(5,3)=EpsC;

% for iu = 1:9
%     UnitCell(iu, 3) = 1 / (iu * 100);
% end


%% Create k points
% k = kPoints(xPeriod, yPeriod, DPoints); %(Rectangular path)
k = kPointsReduced(P1, P2, DPoints); % Square path in k-space
k = cat(1, k(size(k, 1), :), k); % Add gamma point to begining of k-vector
Nk = size(k, 1); % Number of k vectors

% k = [1, 2];
% Nk = 1;

% %plots path
% hold off;
% for i = 1:Nk
%     scatter(k(i, 1), k(i, 2));
%     pause(0.001);
%     hold on;
%     axis([min(k(:, 1)) max(k(:, 1)) min(k(:, 2)) max(k(:, 2))])
%     axis square
% end

%% Create cell containing the A matrices
% FDFDmatrix is an array of Nk number of A matrices (in the generalized eigenvalue problem)
% There is a separate one for each k vector.
% UnitRow = 1 ./ UnitCell(:, 3);
% p = zeros(sqPix); %Setting up a matrix to contain the dielectric values at positions corresponding to their mesh positions in the unit cell
% for i = 1:sqPix
%     p(i, :) = UnitRow(((i - 1) * sqPix + 1):(i * sqPix));
% end
% FDFDmatrix = TEMatrix(P1, P2, xPix, yPix, k, p);
[FDFDmatrix, p] = TEMatrix2(P1, P2, xPix, yPix, k, UnitCell);
% Real = real(FDFDmatrix{1});
% Imag = imag(FDFDmatrix{1});
% FDFDmatrix = TEMatrix(P1, P2, xPix, yPix, k, UnitCell);
% See notes above and https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem

%% Solve for eigenvalues
eValues = cell(Nk, 1); % Initializing matrices
eModes = cell(Nk, 1);
eModesNorm = cell(Nk, 1);
for m = 1:Nk
    eValues{m} = eig(FDFDmatrix{m}); % Calculates Nmesh eigenvalues for each k point
    eModes{m} = sqrt(eValues{m}) * c;  % Calculates the eModes = omega (angular frequency)
    eModesNorm{m} = (eModes{m} * P1) / (2 * pi * c); % Normalizes the frequency
    lbdaModes{m} = (2 * pi * c ./ eModes{m}) * 1e9; % Solve for the wavelengths to use in COMSOL comparisons
end   
eModesNormMatrix = zeros(Nk, T * Nmesh); %Instead of a Nk * 1 * Nmesh (3D) matrix, simplifying to a Nk * Nmesh (2D) matrix to contain normalized freq values
lbda = zeros(Nk, T * Nmesh); % Instead of a Nk * 1 * Nmesh (3D) matrix, making a Nk * Nmesh (2D) matrix to contain wavelength values 
for j = 1:Nk % Cycling through all k points
    eModesNormMatrix(j, :) = eModesNorm{j}'; % Matrix with a row of normalized frequency values for each k point
    lbda(j, :) = lbdaModes{j}'; % Matrix with a row of wavelength values for each k point
end

%% plot band structure
kPlotVector = zeros(Nk, 1); % Initializes x-axis (k)
for i = 1:Nk
    kPlotVector(i) = i;
end
%   eModesNormMatrix = cat(2, kPlotVector, eModesNormMatrix); 
%   PlotModes = plot(eModesNormMatrix(:, 1), eModesNormMatrix(:, 2:(xPix * yPix + 1)));
%   PlotModes = plot(eModesNormMatrix(:, 1), eModesNormMatrix(:, 2:Nmesh));


PlotModes = zeros(Nk, 1);
for i = 1:(T * Nmesh)
% for i = (Nmesh + 1):(2 * Nmesh) % Cycle through all of the solved eigenvalues for a given k
    PlotModes(i) = plot(kPlotVector, eModesNormMatrix(:, i), '.r');
    hold on;
end

% Create x-axis labels for high symmetry k-points
Xpoint = round(DPoints / sqrt(2));
set(gca, 'Xtick', [1, (round(DPoints / sqrt(2)) + 1), (2 * round(DPoints / sqrt(2)) + 1), Nk]); 
set(gca, 'XTickLabel', ['G'; 'X'; 'S'; 'G']);
% %     text(1, -0.075 * 10^16, '\Gamma');
% %     text(Nk, -0.075 * 10^16, '\Gamma');
axis([1 Nk 0 wmax])
hold on 

% Make the vertical lines for high symmetry k-points
plot([(Xpoint + 1) (Xpoint + 1)], [0 wmax], 'k');
plot([(2 * Xpoint + 1) (2 * Xpoint + 1)], [0 wmax], 'k');
hold off

%% Plot band structure and COMSOL results together for comparison
% figure
% % lbda = real((2 * pi * c ./ eModesNormMatrix) * Period); % Solve for the wavelengths to use in COMSOL comparisons
%if EpsB == 1
%    ComsolComparePosts(Nk, Nmesh, lbda, Xpoint);
%elseif EpsC == 1
%    ComsolCompareHoles(Nk, Nmesh, lbda, Xpoint);
%end

%% Plot matrix via spy for testing
% figure
% Pix = sqPix;
% spy(FDFDmatrix{1})
% hold on
% plot([Pix*Pix + 0.5 Pix*Pix + 0.5],[0 Pix*Pix*2+ 0.5],'k-');
% plot([0 Pix*Pix*2+ 0.5],[Pix*Pix + 0.5 Pix*Pix + 0.5],'k-');
% hold off

disp('Code has finished running.');