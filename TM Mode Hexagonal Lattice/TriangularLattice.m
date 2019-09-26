A1 = [100, 0];
A2 = [50, sqrt(3) * 50];

for i = 1:50
    for j = 1:50
        M{i, j} = (i * A1) + (j * A2);
    end
end