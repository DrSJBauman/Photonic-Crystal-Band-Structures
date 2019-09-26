function [k] = kPointsReducedRegHexLattice(Period, DPoints)
% input: Period = length in m (length of rhombus unit cell side). 
% input: Dpoints = length of diagonal (hypotenuse) of irreducible brillouin zone.

%Refer to the powerpoint file in this Hexagonal TM Mode folder to visualize the Brillouin Zone (BZ)

    %% Unit vectors of lattice in real space, a1 and a2, stored as rows of a matrix A
    %A = [Period, Period / 2; 0, sqrt(3) * Period / 2]; %original code
    A = [Period, 0; Period / 2, sqrt(3) * Period / 2]; %2X2 matrix of the basis vectors

    % Unit vectors of lattice in k-space
    B = (2 * pi * eye(2)) / A; %The formula for magnitude of the k-space vectors
    B = transpose(B); %Switching rows and columns
    
    %% First irreducible Brillouin zone points. These were hard coded vs. using B from above
    Gamma = [0; 0]; %Gamma is at the center of the BZ
%    M = ( pi / Period  ) * [sqrt(3) / 2 ; 1/2 ]; %updated from Scheffler et al Ch 5 "Theoretical Material Science"
    M = ( pi / (Period * sqrt(3)) ) * [sqrt(3) ; -1 ]; %updated from Salva
%    K = M + [(B(2, 2) / (sqrt(3) * 2)); 0];
%    K = ( 2 * pi / ( Period ) ) * [ 1/sqrt(3) ; 0 ]; %updated from Scheffler et al 
    K = ( 4 * pi / ( Period * 3) ) * [ 1 ; 0 ]; %updated from Salva 
    
    %DPoints = Number of k points along Gamma to K
    Ngm = sqrt(3) * DPoints / 2; % Number of k points along Gamma to M
    Nmk = DPoints / 2; % Number of k points along M to K
    Ngm = round(Ngm); %Round up to next highest integer
    Nmk = round(Nmk);
    
    %% k points along path from Gamma to M
    k_GM = zeros(Ngm, 2); %Setting up a matrix that is Ngm rows by 2 columns
    for i = 1:Ngm
        %The columns are set depending on this function that incorporates Gamma and M
        k_GM(i, 1) = (Gamma(1, 1) * (1 - (i / Ngm))) + (M(1, 1) * (i / Ngm));
        k_GM(i, 2) = (Gamma(2, 1) * (1 - (i / Ngm))) + (M(2, 1) * (i / Ngm));
    end
    
    % k points along path from M to K 
    k_MK = zeros(Nmk, 2); %Setting up a matrix that is Nmk rows by 2 columns
    for i = 1:Nmk
        %The columns are set depending on this function that incorporates M and K
        k_MK(i, 1) = (M(1, 1) * (1 - (i / Nmk))) + (K(1, 1) * (i / Nmk));
        k_MK(i, 2) = (M(2, 1) * (1 - (i / Nmk))) + (K(2, 1) * (i / Nmk));
    end
    
    % k points along path from K to Gamma
    k_KG = zeros(DPoints, 2); %Setting up a matrix that is DPoints rows by 2 columns
    for i = 1:DPoints
        %The columns are set depending on this function that incorporates K and Gamma
        k_KG(i, 1) = (K(1, 1) * (1 - (i / DPoints))) + (Gamma(1, 1) * (i / DPoints));
        k_KG(i, 2) = (K(2, 1) * (1 - (i / DPoints))) + (Gamma(2, 1) * (i / DPoints));
    end
    k = [k_GM; k_MK; k_KG]; %Putting the k vectors into one matrix
end


