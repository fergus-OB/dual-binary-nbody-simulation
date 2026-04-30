function F = MyDE(~, yy)
    % Define the masses of the stars (m1 and m2)
    m1 = 0.9;  % mass of star 1 (arbitrary units)
    m2 = 1.1;  % mass of star 2 (arbitrary units)
    
    % Define the gravitational constant (G)
    G = 1;  % normalized gravitational constant (arbitrary units)
    
    % Unpack the state vector yy into the positions and velocities of the stars
    x1 = yy(1);   % position of star 1 in x
    y1 = yy(2);   % position of star 1 in y
    vx1 = yy(3);  % velocity of star 1 in x
    vy1 = yy(4);  % velocity of star 1 in y
    
    x2 = yy(5);   % position of star 2 in x
    y2 = yy(6);   % position of star 2 in y
    vx2 = yy(7);  % velocity of star 2 in x
    vy2 = yy(8);  % velocity of star 2 in y
    
    % Calculate the distance between the two stars
    r12 = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    
    % Equations of motion for the binary system (Newton's law of gravitation)
    ax1 = -G * m2 * (x1 - x2) / r12^3;  % acceleration of star 1 in x
    ay1 = -G * m2 * (y1 - y2) / r12^3;  % acceleration of star 1 in y
    
    ax2 = -G * m1 * (x2 - x1) / r12^3;  % acceleration of star 2 in x
    ay2 = -G * m1 * (y2 - y1) / r12^3;  % acceleration of star 2 in y

    % Kinetic energy
    KE1 = 0.5 * m1 * (vx1^2 + vy1^2);  % for star 1
    KE2 = 0.5 * m2 * (vx2^2 + vy2^2);  % for star 2
    KE_total = KE1 + KE2;

    % Potential energy (gravitational potential)
    PE = -G * m1 * m2 / r12;

    % Total energy
    E_total = KE_total + PE;

    %disp(['Total Energy: ', num2str(E_total)]);
    
    % Return the derivatives as a column vector
    F = [vx1; vy1; ax1; ay1; vx2; vy2; ax2; ay2];
end

% Main script
tend = 1000;  % End time for the simulation

% Initial conditions (positions and velocities)
ICs = [1, 0, 0, 0.5, -1, 0, 0, -0.5];  % example for stable orbital motion
  % Example initial conditions

RTOL = 1E-10; ATOL = 1E-10;
options = odeset('RelTol', RTOL, 'AbsTol', ATOL);

% Solve the equations using ode45
[t, YE] = ode45(@MyDE, [0, tend], ICs, options);

% Calculate the relative distance between the two stars
distance = sqrt((YE(:,5) - YE(:,1)).^2 + (YE(:,6) - YE(:,2)).^2);

% Plot the relative distance between the stars over time
figure;
plot(t, distance);
xlabel('Time');
ylabel('Relative Distance');
title('Relative Distance Between Two Stars Over Time');

% Optional: Phase plot of the positions of the stars
figure;
plot(YE(:,1), YE(:,2), 'r', 'DisplayName', 'Star 1 Orbit');
hold on;
plot(YE(:,5), YE(:,6), 'b', 'DisplayName', 'Star 2 Orbit');
xlabel('X Position');
ylabel('Y Position');
title('Phase Plot of Star Orbits');
%disp(['Total Energy: ', num2str(E_total)]);
legend;
