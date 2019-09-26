function [FDFDmatrix] = kMatrixRegHex(Period, Pix, k)
%input Period: Period of the photonic crystal
%input Pix: Total Number of mesh elements in one direction
%input k: The k-path values received from the kPoints function
%output FDFDmatrix is an array of A matrices, one for each k point. note

%% Building the k-matrix that is contained in the PPt/pdf in this Hexagonal TM Mode folder
    B = [Period, 0; Period / 2, sqrt(3) * Period / 2]; %basis vectors (A) of triangular lattice
    a = Period / Pix; % size of one mesh element (a=dx=dy)
    FDFDmatrix = cell(size(k, 1), 1); % Initialize the function's output matrix
    
    for m = 1:size(k, 1) % iterates through each k point
        KA1{m} = k(m, 1) * B(1, 1) + k(m, 2) * B(1, 2); %k dot A1
        KA2{m} = k(m, 1) * B(2, 1) + k(m, 2) * B(2, 2); %k dot A2
        %% Main diagonal
        Main_Diag{m} = zeros(Pix); %Initializing the Pix-by-Pix matrix that repeats along the main diagonal of the large matrix
        for i = 1:Pix
            Main_Diag{m}(i, i) = 6 / a^2; %Diagonal positions all the way down (1,1) (2,2) etc
        end
        for i = 2:Pix
            Main_Diag{m}(i - 1, i) = - 1 / a^2; %First row, second column etc. (1,2) (2,3) etc
            Main_Diag{m}(i, i - 1) = - 1 / a^2; %Second row, first column etc. (2,1) (3,2) ... 
        end
        % next line: for (1,3) for 3x3 mesh
        Main_Diag{m}(1, Pix) = - (1 / (a^2)) * exp(1i * KA2{m}); %First row, last column in the Pix-by-Pix matrix (with paranthesis)
        % next line: for (3,1) for 3x3 mesh
        Main_Diag{m}(Pix, 1) = conj(Main_Diag{m}(1, Pix)); %Last row, first column
        
        %% Right diagonal
        Right_Diag{m} = zeros(Pix);
        for i = 1:Pix %fills in the diagonal
            Right_Diag{m}(i, i) = - 1 / a^2; %3x3: (1,4) (2,5) ...
        end
        for i = 2:Pix % second row first col, etc. 3x3: (2,4)
            Right_Diag{m}(i, i - 1) = - 1 / a^2;
        end
        %next line: for 3x3: (1,6)
         Right_Diag{m}(1, Pix) = - (1 / (a^2)) * exp(1i * KA2{m}); 

        %% Left diagonal
        for i = 1:Pix %3x3: (4,1) (5,2) ...
            Left_Diag{m}(i, i) = - 1 / a^2;
        end
        for i = 2:Pix %3x3: (4,2) (5,3) ...
            Left_Diag{m}(i - 1, i) = - 1 / a^2;
        end
        Left_Diag{m}(Pix, 1) = conj(Right_Diag{m}(1, Pix)); %3x3: (6,1)
        
        %% Top diagonal (Top Right Matrix)
        Top_Diag{m} = zeros(Pix);
        for i = 1:Pix % 3x3: (1,7) (2,8) ...
            Top_Diag{m}(i, i) = - (1 / a^2) * exp(1i * KA1{m});
        end
        for i = 2:Pix % 3x3: (1,8) (2,9) ...
            Top_Diag{m}(i - 1, i) = - ((1 / a^2) * exp(1i * KA1{m}));
        end
        Top_Diag{m}(Pix, 1) = - (1 / a^2) * exp(1i * (- KA1{m} + KA2{m}));
        
        %% Bottom diagonal (Bottom Left Matrix)
        Bot_Diag{m} = zeros(Pix);
        for i = 1:Pix %Diagonal 3x3: (7,1) (8,2) ...
            Bot_Diag{m}(i, i) = conj(Top_Diag{m}(i, i));
        end
        for i = 2:Pix %left diag 3x3: (8,1) (9,2) ...
            Bot_Diag{m}(i, i - 1) = conj(Top_Diag{m}(i - 1, i));
        end 
        Bot_Diag{m}(1, Pix) = conj(Top_Diag{m}(Pix, 1)); %3x3: (7,3)
        
        %% Create empty unit cell matrix
        zeroUnitCell = zeros(Pix);
        % Initialize final matrix
        for j = 1:Pix %initilizes first A to all zeros, each cell is a matrix for now
            for g = 1:Pix
                FDFDmatrix{m}{j, g} = zeroUnitCell;
            end
        end
        % Create cell diagonals
        for j = 1:Pix
            FDFDmatrix{m}{j, j} = Main_Diag{m};
        end
        % Create adjacent cell diagonals
        for j = 1:(Pix - 1)
            FDFDmatrix{m}{j, j + 1} = Right_Diag{m};
        end
        for j = 2:Pix
            FDFDmatrix{m}{j, j - 1} = Left_Diag{m};
        end
        % Create outer matrix
        FDFDmatrix{m}{1, Pix} = Top_Diag{m};
        FDFDmatrix{m}{Pix, 1} = Bot_Diag{m};
        % Convert final matrix cell array to array
        FDFDmatrix{m} = cell2mat(FDFDmatrix{m});
    end
end

