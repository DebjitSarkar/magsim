%% MAGSIM
% Debjit Sarkar
% The purpose of this script is to simulate the electromagnetic forces
% between point charges as an intermediate step to doing FEM/FDTD analysis of
% larger systems like transistors

% Note: Does not include second-order field effects yet

%% PHYSICAL CONSTANTS
epsilon0 = 8.854e-12; % permittivity of free space
mu0 = 4 * pi * 1e-7;  % permeability of free space
me = 9.109e-31;       % mass of electron
mp = 1.673e-27;       % mass of proton
e = 1.602e-19;        % fundamental charge
c = 3e8;              % light speed in vacuum

%% SIMULATION PARAMETERS AND VARIABLES
numPoints = 3;        % number of point charges
tbegin = 0;           % simulation start time
step = 0.001;         % time step
numSteps = 10;        % number of time steps
tend = tbegin + numSteps * step; % simulation end time

debugLvl = 1; % 0 = off, 1 = time loop, 2 = calc loop

%% POINT CHARGES
charges = zeros(numPoints, 13); 
%   Field:  ID, charge, mass, position, velocity, acceleration, mobility
%   Index:  1   2       3     4  5  6   7  8  9   10 11 12      13
% Example: [1,  1e-10,  1e-9, 0, 0, 0,  5, 0, 0,  0, 0, -9.8,   1]
% note about mobility: 1 = mobile, 0 = static
charges(1, :) = [1, -e, me, 0 0 0, 0 0 0, 0 0 0, 1]; % electron
charges(2, :) = [2, -e, me, 0 1 0, 0 0 0, 0 0 0, 1]; % electron
charges(3, :) = [3,  e, mp, 1 0 0, 0 0 0, 0 0 0, 0]; % proton

%% TEST CASES
% Remember to change the number of points

% Test for electric attraction and repulsion
%charges(1, :) = [1, -e, me, 0 0 0, 0 0 0, 0 0 0]; % electron
%charges(2, :) = [2, -e, me, 0 1 0, 0 0 0, 0 0 0]; % electron
%charges(3, :) = [3,  e, mp, 1 0 0, 0 0 0, 0 0 0]; % proton

% Test for magnetic attraction and repulsion
%charges(1, :) = [1, -e, me, 0 0 0,   0 1e10 0, 0 0 0]; % electron
%charges(2, :) = [2, -e, me, 100 0 0, 0 1e10 0, 0 0 0]; % electron
%charges(3, :) = [3,  e, mp, 50 0 0,  0 1e10 0, 0 0 0]; % proton

%% CALCULATIONS
figure;
hold on;

for t = tbegin:step:tend
    charges(:, 10:12) = zeros(size(charges, 1), 3);
    if debugLvl >= 1
        disp("Velocities at t = " + t);
        disp(charges(:, 7:9));
    end
    for i = 1:numPoints
        if charges(i, 13) == 1 % Skips calculations for static charges
            for j = 1:numPoints
                if(i ~= j)
                    if debugLvl >= 2
                        disp("t = " + t + "|i = " + i + "|j = " + j);
                    end
                    posDiff = charges(j, 4:6) - charges(i, 4:6);
                    posNorm = norm(posDiff);
                    %F = q * (E + v x B)
                    %a = (q / m) * (E + v x B)
                    %B = (mu0 / (4 * pi)) * (q / r^3) * cross(v, r)
                    B = (mu0 * charges(j, 2) * ...
                        cross(charges(j, 7:9), posDiff)) / (4 * pi * posNorm^3);
                    %B = [0 0 0];
                    E = (charges(j, 2) * posDiff) / ...
                        (4 * pi * epsilon0 * posNorm^3);
                    %a = (qE + qv x B) / m
                    % TODO: verify that the acceleration is in the right
                    % direction (the negative signs in front of E and
                    % cross(...)
                    charges(i, 10:12) = charges(i, 10:12) + (charges(i, 2) /...
                        charges(i, 3)) * (-E + -cross(charges(i, 7:9), B));
                end
            end
        end
    end
    for i = 1:numPoints
        if charges(i, 13) == 1
            charges(i, 7:9) = charges(i, 7:9) + step * charges(i, 10:12);
            charges(i, 4:6) = charges(i, 4:6) + step * charges(i, 7:9);
        end
    end
    
%% PLOTTING
    % plots charges from white = beginning to black = end
    colorSpec = [1 1 1] * (1 - (t - tbegin) / tend);
    plot3(charges(:, 4), charges(:, 5), charges(:, 6), 'o', ...
        'MarkerFaceColor', colorSpec, ...
        'MarkerEdgeColor', 'k'); % k = black
end
