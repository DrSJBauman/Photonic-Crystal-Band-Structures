function [k] = kPointsReducedRegHexLattice(Period, DPoints)
    % Unit vectors of lattice in real space
    A = [Period, Period / 2; 0, sqrt(3) * Period / 2]; %original code
    %A = [Period, 0; Period / 2, sqrt(3) * Period / 2]; %edited for test
    % Unit vectors of lattice in k-space
    B = (2 * pi * eye(2)) / A;
    B = transpose(B);
    % First irreducible Brillouin zone points
    Gamma = [0; 0];
    M = B(:, 2) / 2;
    K = M + [(B(2, 2) / (sqrt(3) * 2)); 0];
    % Number of k points along the horizontal (Gamma to X) and vertical (Gamma to y)
    Nx = sqrt(DPoints^2 / (1 + (B(2, 2) / B(1, 1))^2));
    Ny = sqrt(DPoints^2 - Nx^2);
    Nx = round(Nx);
    Ny = round(Ny);
    % % k points along path from Gamma to M
    k_GM = zeros(Ny, 2);
    for i = 1:Ny
        k_GM(i, 1) = (Gamma(1, 1) * (1 - (i / Ny))) + (M(1, 1) * (i / Ny));
        k_GM(i, 2) = (Gamma(2, 1) * (1 - (i / Ny))) + (M(2, 1) * (i / Ny));
    end
    % % k points along path from M to K 
    k_MK = zeros(Nx, 2);
    for i = 1:Nx
        k_MK(i, 1) = (M(1, 1) * (1 - (i / Nx))) + (K(1, 1) * (i / Nx));
        k_MK(i, 2) = (M(2, 1) * (1 - (i / Nx))) + (K(2, 1) * (i / Nx));
    end
    % % k points along path from K to Gamma
    k_KG = zeros(DPoints, 2);
    for i = 1:DPoints
        k_KG(i, 1) = (K(1, 1) * (1 - (i / DPoints))) + (Gamma(1, 1) * (i / DPoints));
        k_KG(i, 2) = (K(2, 1) * (1 - (i / DPoints))) + (Gamma(2, 1) * (i / DPoints));
    end
    k = [k_GM; k_MK; k_KG];
end


