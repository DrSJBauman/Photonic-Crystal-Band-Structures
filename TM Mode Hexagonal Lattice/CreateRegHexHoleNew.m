function [HexMatrix, radius] = CreateRegHexHoleNew(xPix, DielectConst, rRatio)
%Creates a single unit cell for a hexagonal lattice with air holes
%surrounded by a material with the dielectric constant provided.
    
Pix=xPix;
eps = DielectConst;

%xPix: number of mesh elements 
radius = rRatio * xPix; %radius
r = radius;

rhom1 = zeros(Pix,round(1.5 * Pix)); %a rectangular matrix which can contain the rhombus unit cell
sqPx = eps*ones(Pix,Pix);

L = 1.5*Pix;

for i=1:Pix %loop for row
    for j=1:Pix %loop for each element in row
        rhom1(i,round(j+0.5*(Pix-i+1)))=sqPx(i,j);
    end
end

%% Draw circles

%TL circle

cx=Pix/2;
cy=0;

for i=1:Pix %loop for row
    for j=1:L %loop for each element in row
        if i<( sqrt(r^2 - (j-cx)^2) +cy)
            rhom1(i,j)=1;
        end
    end
end

%TR circle

cx=L;
cy=0;

for i=1:Pix %loop for row
    for j=1:L %loop for each element in row
        if i< sqrt(r^2 - (j-cx)^2) +cy
            rhom1(i,j)=1;
        end
    end
end

cx=0;
cy=Pix;

for i=1:Pix %loop for row
    for j=1:L %loop for each element in row
        if i> ( - sqrt(r^2 - (j-cx)^2) + cy )
            rhom1(i,j)=1;
        end
    end
end


cx=Pix;
cy=Pix;

for i=1:Pix %loop for row
    for j=1:L %loop for each element in row
        if i> ( - sqrt(r^2 - (j-cx)^2) + cy )
            rhom1(i,j)=1;
        end
    end
end

 %% Shift back to square matrix 

HexMatrix = zeros(Pix,Pix);
 
for i=1:Pix %loop for row
    for j=1:Pix %loop for each element in row
        HexMatrix(i,j)=rhom1(i,round(j+0.5*(Pix-i+1)));
    end
end

HexMatrix = flipud(HexMatrix);

%  %ax=pcolor(rhom1);
%  ax=pcolor(HexMatrix);
%  load sjmap.mat
%  colormap(sjmap);
%  set(ax, 'EdgeColor', 'none');
 
end