%% This function creates the unit cell for any 2D geometry by utilizing the two basis vectors.
% It takes the following inputs in order to create a (mesh size squared) by 3 matrix 
% containing the coordinates and dielectric constant values for each point in the mesh.

% This resolves the issue of the previous method in which we had tried to
% use a matrix as the unit cell itself, limiting us to a square structure
% that had been sheared from a diagonal lattice.

% In this new geometry, the mesh will actually have unused coordinates in
% two of the corners.

clear all
clc
%% Inputs
P1 = 100 * 10^-9; %Period along a1 [m]
P2 = 100 * 10^-9; %Period along a2 [m]
theta = 30; %Angle between a1 and a2 vectors from 0 to 90 [degrees]
n1 = 21; %Number of mesh elements along a1
n2 = 21; %Number of mesh elements along a2
Eps1 = 13; %Dielectric constant of posts/holes [units]
Eps2 = 1; %Dielectric constant of bulk [units]
rRatio = 0.25; %Radius divided by the center-to-center distance between posts/holes (choose P1 or P2)

%% Setting things up
r = rRatio * P1; %Radius of circular post/hole [m]
a1 = [1, 0]; %Basis vector 1 (multiply by P1)
a2 = [cos(theta / 180 * pi), sin(theta / 180 * pi)]; %Basis vector 2 (multiply by P2)
L = ceil(n1 + (n2 / tan(theta * pi/180))); %Total mesh elements along a1 in the rectangular matrix
Geometry = zeros(n2, L); %Rectangular matrix which can contain the rhombus unit cell
RhomValues = Eps2*ones(n2,n1); %Rectangular matrix containing Eps2 values at all coordinates
UnitCell = zeros(n1 * n2, 3); %Matrix that will contain the coordinates and Eps values for each mesh point

%% Creating the Unit Cell

for i=1:n2 %Loop for row
    for j=1:n1 %Loop for each column in RhomValues
%         Geometry(i * round(a2(2)),j + i * round(a2(1))) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions
%         Geometry(i, round(j - n2 + 8 + (L - i))) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions
%         Geometry(i, round(i * tan(theta))+1) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions
%         Geometry(i, ceil(j + (tan(theta) * (L-i)))) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions
        Geometry(i, (ceil(tan(theta / 180 * pi) * i)+ (L - j))) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions

%         for k=1:L %Loop for each column in rectangular matrix Geometry
%             if (L - n1 - i - 1) <= k 
%                 if k <= (L - i - 1) %Setting the borders of the rhombus
%                     Geometry(i,k) = RhomValues(i,j); %Assigns Eps2 values to rhombus positions
%                 end
%             else Geometry(i,k) = 0;    
%             end
%         end
    end
end
% 
% for i=1:n2 %loop for each row
%     for j=1:n1 %loop for each element in row
%         if UnitCell(i,j in BotLeftCircle
%             UnitCell(i,j) = [a1 * i, a2 * j, Eps1]
%     end
% end
%         
% % Top Left Circle
% Ci = a1*(n1 / 2); %Coordinate of the circle center in the a1 direction
% Cj = ;
%         
% 
% %% Output
% UnitCell = 
