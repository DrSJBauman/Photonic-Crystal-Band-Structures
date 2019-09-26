function [UnitCell] = GenGeom(P1, P2, theta, n1, n2, EpsC, EpsB, rRatio)
%% This function creates the unit cell for any 2D geometry by utilizing the two basis vectors.
% It takes the following inputs in order to create a (mesh size squared) by 3 matrix 
% containing the coordinates and dielectric constant values for each point in the mesh.

% This resolves the issue of the previous method in which we had tried to
% use a matrix as the unit cell itself, limiting us to a square structure
% that had been sheared from a diagonal lattice.

%% Inputs
%P1 = Period along a1 [nm]
%P2 = Period along a2 [nm]
%theta = Angle between a1 and a2 vectors from 0 to 90 [degrees]
radtheta = theta * pi / 180; %theta converted to [radians]
%n1 = Number of mesh elements along a1 [elements]
%n2 = Number of mesh elements along a2 [elements]
%EpsC = Dielectric constant of the cylinders (posts/holes) [unitless]
%EpsB = Dielectric constant of the bulk [unitless]
%rRatio = Radius divided by the center-to-center distance between posts/holes (choose P1 or P2)

%% Setting things up
r = rRatio * P1; %Radius of circular post/hole [m]
a1 = P1/n1 * [1, 0]; %Basis vector 1 [nm/mesh]
a2 = P2/n2 * [cos(radtheta), sin(radtheta)]; %Basis vector 2 [nm/mesh]
UnitCell = zeros(n1 * n2, 3); %Matrix that will contain the coordinates and Eps values for each mesh point

%% Creating the matrix containing the Unit Cell coordinates and dielectric values
ni = 1; %Initialize counter
for i=1:n1 %Loop for row in the n1*n2 space
    for j=1:n2 %Loop for each column in the n1*n2 space
         x = (j - 1) * a1(1) + (i - 1) * a2(1); %Establish x-coordinates using basis vector dot product
         y = (j - 1) * a1(2) + (i - 1) * a2(2); %Establish y-coordinates using basis vector dot product
         
         if (x - cos(radtheta) * P1)^2 + (y - sin(radtheta)*P2)^2 < r^2 %inside TL circle
            UnitCell(ni, :) = [x y EpsC];
            color='r';
         elseif (x - (P1 + cos(radtheta) * P1))^2 + (y - sin(radtheta)*P2)^2 < r^2 %inside TR circle
            UnitCell(ni, :) = [x y EpsC];
            color='r';
         elseif x^2 + y^2 < r^2 %inside BL circle
            UnitCell(ni, :) = [x y EpsC];
            color='r';
         elseif (x-P1)^2 + y^2 < r^2 %inside BR circle
            UnitCell(ni, :) = [x y EpsC];
            color='r';
         else 
            UnitCell(ni, :) = [x y EpsB];
            color='b';
          end
         ni = ni + 1;

%         Plotting the values of UnitCell, one at a time, to visualize the unit cell
        plot(x,y,'s','MarkerSize',3,'MarkerEdgeColor',color,'MarkerFaceColor',color);
        hold on
    end
end
% %% Making the plot look pretty
axis ([0 max(max(x), max(y)) 0 max(max(x), max(y))]);
axis square;

end

