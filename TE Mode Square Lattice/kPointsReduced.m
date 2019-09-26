function [k] = kPointsReduced(xPeriod, yPeriod, DPoints)
    % Unit vectors of lattie in real space
    A = [xPeriod, 0; 0, yPeriod];
    % Unit vectors of lattic in k-space
    B = (2 * pi * eye(2)) / A;
    B = transpose(B);
    % First irreducible Brillouin zone points
    Gamma = [0; 0];
    S = (B(:, 1) + B(:, 2)) / 2;
    X = B(:, 1) / 2;
    % Number of k points along the horizontal (Gamma to X) and vertical (Gamma to y)
    Nx = sqrt(DPoints^2 / (1 + (B(2, 2) / B(1, 1))^2));
    Ny = sqrt(DPoints^2 - Nx^2);
    Nx = round(Nx);
    Ny = round(Ny);
    % k points along path from Gamma to X
    k_GX = zeros(Nx, 2);
    for i = 1:Nx
        k_GX(i, 1) = (Gamma(1, 1) * (1 - (i / Nx))) + (X(1, 1) * (i / Nx));
        k_GX(i, 2) = (Gamma(2, 1) * (1 - (i / Nx))) + (X(2, 1) * (i / Nx));
    end
    % k points along path from X to S (S = M)
    k_XS = zeros(Ny, 2);
    for i = 1:Ny
        k_XS(i, 1) = (X(1, 1) * (1 - (i / Ny))) + (S(1, 1) * (i / Ny));
        k_XS(i, 2) = (X(2, 1) * (1 - (i / Ny))) + (S(2, 1) * (i / Ny));
    end
    % k points along path from Gamma to S
    k_SG = zeros(DPoints, 2);
    for i = 1:DPoints
        k_SG(i, 1) = (S(1, 1) * (1 - (i / DPoints))) + (Gamma(1, 1) * (i / DPoints));
        k_SG(i, 2) = (S(2, 1) * (1 - (i / DPoints))) + (Gamma(2, 1) * (i / DPoints));
    end
    k = [k_GX; k_XS; k_SG];
end


