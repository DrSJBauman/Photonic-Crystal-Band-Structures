function lbdaVSpix = GapShiftTest(Nk,Pix, Nmesh, lbda, Mpoint, Kpoint)
%GapShiftTest plots the two closest wavelength values of the main gap from
%the FDFD code versus number of mesh elements so we can see if the number 
%of mesh elements is causing a shift in the results.

%% Inputs
%Nk = Number of total k points taken from the band structure main code
%Pix = Number of mesh elements along one direction
%Nmesh = Number of total mesh points from the band structure main code(Pix^2)
%lbda = (Nk*Nmesh) Matrix of wavelengths from the eigenstates calculated in main code (nm)
%MPoint = k-coordinates of the high symmetry point in the irreducible BZ
%KPoint = k-coordinates of the high symmetry point in the irreducible BZ

%% Capturing the desired wavelength values for the lowest energy gap at M
RowM = zeros(1, Nmesh); %Initializing vector M
RowM = lbda(Mpoint, :); %Selecting all wavelengths along M
GapPoints = cat(1, importdata('GapPoints.mat'), [RowM(1), RowM(2), Pix]); %Storing the gap points in a matrix for each different mesh size 
save('GapPoints.mat', 'GapPoints'); %Saving the result so that it isn't cleared during each new run with different mesh sizes

%% Plotting wavelength vs # of mesh elements for the lowest energy gap at M
for i = 1:Pix %Count up to the number of mesh elements
    lbdaVSpix = plot(GapPoints(:,3), GapPoints(:,1), 'bo', GapPoints(:,3), GapPoints(:,2), 'go');
    title('Wavelength shift for different numbers of mesh elements');
    ylabel('Wavelength (nm)');
    xlabel('Mesh elements');
    legend('1st band','2nd band')
end

