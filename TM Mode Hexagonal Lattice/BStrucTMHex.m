clear all
clc
clf

%% Parameters:    
Period = 200 * 10^-9; % Distance between atoms [m]
P1 = Period; %Distance along basis vector 1 [m]
P2 = Period; %Distance along basis vector 2 [m]
DPoints = 25; % Number of k points along the diagonal of the calculated irreducible BZ
Pix = 21; % Number of mesh points in the photonic crystal unit cell in one direction
%NOTE: Fix geometry - dielectric smoothing?
Pix = Pix; % Sloppy coding?
wmax = 0.8; % Maximum y-axis value on the final band structure plot
n1 = Pix; % Number of mesh points along basis vector 1
n2 = Pix; % Number of mesh points along basis vector 2
theta = 60; % Angle between the real space basis vectors [degrees]
%GaAs has dielectric constant of ~11.8?
EpsC = 11.8; % Relative dielectric constant of the cylinders (posts/holes) [unitless]
EpsB = 1; % Relative dielectric constant of the bulk surrounding the posts/holes [unitless]
rRatio = 0.25; % Ratio between the radius of the cylinders and the period 
if rRatio >= 0.5
    error('Try smaller ratio');
end
c =  2.99792458 * 10^8; % Speed of light [m/s]
Nmesh = Pix^2; % Number of total mesh elements in the unit cell in real space

%% Create unit cell by applying dielectric constant values to mesh positions.
UnitCell = GenGeom(P1, P2, theta, n1, n2, EpsC, EpsB, rRatio);

%% Create epsilon matrix
epsA = zeros(Pix^2);
for i = 1:Pix^2
    epsA(i, i) = UnitCell(i, 3); %Creates a diagonal matrix of epsilon values.
end

%% Create k points
k = kPointsReducedRegHexLatticeJBH(Period, DPoints); % Hex path in k-space
k = cat(1, k(size(k, 1), :), k); % Add gamma point to begining of k-vector
Nk = size(k, 1); % Number of k vectors

%% plots k - path
% hold off;
% for i = 1:Nk
%     scatter(k(i, 1), k(i, 2));
%     pause(0.001);
%     hold on;
%     axis([min(k(:, 1)) max(k(:, 1)) min(k(:, 1)) max(k(:, 1))])
%     axis square
% end

%% Create cell containing the A matrices
% FDFDmatrix is an array of Nk number of A matrices (in the generalized eigenvalue problem)
% There is a separate one for each k vector.
FDFDmatrix = kMatrixRegHex(Period, Pix, k);
% See notes above and https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem

%% Solve for eigenvalues
eValues = cell(Nk, 1); % Initializing matrices
eModes = cell(Nk, 1);
eModesNorm = cell(Nk, 1);
for m = 1:Nk
    eValues{m} = eig(FDFDmatrix{m}, epsA); % Calculates Nmesh eigenvalues = (omega/c)^2 for each k point
    eModes{m} = sqrt(eValues{m}) * c;  % Calculates the eModes = omega (angular frequency)
    eModesNorm{m} = (eModes{m} * Period) / (2 * pi * c); % Normalizes the frequency
%     eModesNorm{m} = eModesNorm{m} * (6 / 7); % Values are shifted by this amount compared to the MIT reference for some reason
% NOTE: Multiplying lbdaModes by (7 / 6) makes the FDFD results align perfectly with COMSOL, and without it, they are shifted a bit.
    lbdaModes{m} = (2 * pi * c ./ eModes{m}) *  1e9; % Solve for the wavelengths to use in COMSOL comparisons
end   
eModesNormMatrix = zeros(Nk, Nmesh); %Instead of a Nk * 1 * Nmesh (3D) matrix, simplifying to a Nk * Nmesh (2D) matrix to contain normalized freq values
lbda = zeros(Nk, Nmesh); % Instead of a Nk * 1 * Nmesh (3D) matrix, making a Nk * Nmesh (2D) matrix to contain wavelength values
for j = 1:Nk % Cycling through all k points
    eModesNormMatrix(j, :) = eModesNorm{j}'; % Matrix with a row of normalized frequency values for each k point
    lbda(j, :) = lbdaModes{j}'; % Matrix with a row of wavelength values for each k point
end

%% plot band structure
kPlotVector = zeros(Nk, 1); % Initializes x-axis (k)
for i = 1:Nk
    kPlotVector(i) = i;
end

PlotModes = zeros(Nk, 1);
for i = 1:(Nmesh) % Cycle through all of the solved eigenvalues for a given k
    PlotModes(i) = plot(kPlotVector, eModesNormMatrix(:, i),'b');
    hold on;
end

% Create x-axis labels for high symmetry k-points
Mpoint = 1 + round( sqrt(3) * DPoints / 2 );
Kpoint = round ( Mpoint + DPoints / 2 );
set(gca, 'Xtick', [1, Mpoint, Kpoint, Nk]); 
set(gca, 'XTickLabel', ['G'; 'M'; 'K'; 'G']);

axis([1 Nk 0 wmax])
hold on 

% % Make the vertical lines for high symmetry k-points
plot([Mpoint Mpoint], [0 wmax], 'k');
plot([Kpoint Kpoint], [0 wmax], 'k');
hold off

%% Plot band structure and COMSOL results together for comparison
figure
% lbda = real((2 * pi * c ./ eModesNormMatrix) * Period); % Solve for the wavelengths to use in COMSOL comparisons
if EpsC == 1
    ComsolCompareHoles(Nk, Nmesh, lbda, Mpoint, Kpoint);
elseif EpsB == 1
    ComsolComparePosts(Nk, Nmesh, lbda, Mpoint, Kpoint);
end

%% Plot to test the shift in wavelength versus mesh size
% figure
% GapShiftTest(Nk, Pix, Nmesh, lbda, Mpoint, Kpoint);

disp('Code has finished running.');