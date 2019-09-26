%% This function creates the unit cell for any 2D geometry by utilizing the two basis vectors.
% It takes the following inputs in order to create a (mesh size squared) by 3 matrix 
% containing the coordinates and dielectric constant values for each point in the mesh.

% This resolves the issue of the previous method in which we had tried to
% use a matrix as the unit cell itself, limiting us to a square structure
% that had been sheared from a diagonal lattice.

% In this new geometry, the mesh will actually have unused coordinates in
% two of the corners.

%% Inputs
P1 = 100; %Period along a1 [nm]
P2 = 100; %Period along a2 [nm]
theta = 60; %Angle between a1 and a2 vectors from 0 to 90 [degrees]
radtheta = theta * pi / 180;
n1 = 40; %Number of mesh elements along a1
n2 = 40; %Number of mesh elements along a2
Eps1 = 13; %Dielectric constant of posts/holes [units]
Eps2 = 1; %Dielectric constant of bulk [units]
rRatio = 0.40; %Radius divided by the center-to-center distance between posts/holes (choose P1 or P2)

%% Setting things up
r = rRatio * P1; %Radius of circular post/hole [m]
a1 = P1/n1 * [1, 0]; %Basis vector 1 (multiply by P1)
a2 = P2/n2 * [cos(radtheta), sin(radtheta)]; %Basis vector 2 (multiply by P2)
L = ceil(n1 + (n2 / tan(radtheta))); %Total mesh elements along a1 in the rectangular matrix
Geometry = zeros(n2, L); %Rectangular matrix which can contain the rhombus unit cell
RhomValues = Eps2*ones(n2,n1); %Rectangular matrix containing Eps2 values at all coordinates
UnitCell = zeros(n1 * n2, 3); %Matrix that will contain the coordinates and Eps values for each mesh point

%% Creating the Unit Cell

ni = 1; %point counter
nc = 1;
for i=1:n1 %Loop for row
    for j=1:n2 %Loop for each column in RhomValues
         %x = UnitCell(ni,1)+a1(1);
         %y = UnitCell(ni,1)+a1(2);
         x = (j - 1) * a1(1) + (i - 1) * a2(1);
         y = (j - 1) * a1(2) + (i - 1) * a2(2);
         if (x - cos(radtheta) * P1)^2 + (y - sin(radtheta)*P2)^2 < r^2 %inside TL circle
            UnitCell(ni, :) = [x y Eps2];
%             color='r';
         elseif (x - (P1 + cos(radtheta) * P1))^2 + (y - sin(radtheta)*P2)^2 < r^2 %inside TR circle
            UnitCell(ni, :) = [x y Eps2];
%             color='r';
         elseif x^2 + y^2 < r^2 %inside BL circle
            UnitCell(ni, :) = [x y Eps2];
%             color='r';
         elseif (x-P1)^2 + y^2 < r^2 %inside BR circle
            UnitCell(ni, :) = [x y Eps2];
%             color='r';
         else 
              UnitCell(ni, :) = [x y Eps1];
%               color='b';
          end
         ni = ni + 1;
         %plot(x,y,color,'MarkerSize',1200000)
%          plot(x,y,'s',...
%             'MarkerSize',3,...
%             'MarkerEdgeColor',color,...
%             'MarkerFaceColor',color);
%          hold on
    end
end

% axis ([0 max(max(x), max(y)) 0 max(max(x), max(y))]);
% axis square;
