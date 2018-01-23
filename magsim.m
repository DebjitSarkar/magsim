%FEM MagSim

epsilon0 = 8.854e-12;
mu0 = 4 * pi * 1e-7;

numPoints = 2;
step = 0.01;
charges = zeros(numPoints, 12); % ID, charge, mass, position, velocity, acceleration
%                         INDEX:  1   2       3     4         7         10
%                       Example: [1,  1e-10,  1e-9, 0, 0, 0,  5, 0, 0,  9.8, 0, 0]

charges(1, :) = [1, 1e-10, 1e-15, 0 0 0, 0 0 0, 0 0 0];
charges(2, :) = [2, 1e-10, 1e-15, 0 1 0, 0 0 0, 0 0 0];
%charges(3, :) = [3, 1e-10, 1e-15, 1 0 0, 0 0 0, 0 0 0];

figure;
hold on;
for t = 0:step:0.1
    for i = 1:numPoints
        for j = 1:numPoints
            if(i ~= j)
                posDiff = charges(j, 4:6) - charges(i, 4:6);
                posNorm = norm(posDiff);
                %F = q * (E + v x B)
                %a = (q / m) * (E + v x B)
                %B = (mu0 / (4 * pi)) * (q / r^3) * cross(v, r)
                %B = (mu0 * charges(j, 2) * cross(charges(j, 7:9), posDiff)) / (4 * pi * posNorm^3);
                B = [0 0 0];
                E = (charges(j, 2) * posDiff) / (4 * pi * epsilon0 * posNorm^3);
                charges(i, 10:12) = (charges(i, 2) / charges(i, 3)) * (E + cross(charges(i, 7:9), B));
            end
        end
        charges(i, 7:9) = charges(i, 7:9) + step * charges(i, 10:12);
        charges(i, 4:6) = charges(i, 4:6) + step * charges(i, 7:9);
    end
    plot3(charges(:, 4), charges(:, 5), charges(:, 6), 'or');
end

%scatter3(charges(:, 4), charges(:, 5), charges(:, 6));
%plot3(charges(1, 4), charges(1, 5), charges(1, 6), 'or');
