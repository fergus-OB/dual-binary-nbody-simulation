function fig4_5
%---------------------------------------------------------
%   Homoclinic-like outer encounter in the HD 74438 2+2 system
%   – shoots for the separatrix that divides “left” vs. “right”
%     passages of the two inner binaries.
%---------------------------------------------------------

%% 0.  Load physical parameters and initial state (24-vector y0)
hd74438_setup          % must define G, m, y0, Pab  (your old block)

%% 1.  Integration options with closest-approach event
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', ...
                 @(t,y)CloseApproach(t,y,G,m));

tgrid   = linspace(0, 80, 501);          % 80 yr shoot window

%% 2.  Two initial conditions that straddle the outer separatrix
ICpos = y0 + [1e-3 0 0  -1e-3 0 0  zeros(1,18)]';  % AB shifted +x
ICneg = y0 - [1e-3 0 0  -1e-3 0 0  zeros(1,18)]';  % AB shifted –x

%% 3.  Interval-bisection shoot
for k = 1:60
    ICmid = 0.5*(ICpos + ICneg);

    [~,~,~,YE] = ode45(@(t,y)QuadrupleRHS(t,y,G,m), ...
                       tgrid, ICmid, options);

    % sign of x-coordinate of AB–CM at first encounter
    rABcm = 0.5*(YE(1,1:3) + YE(1,4:6));
    if rABcm(1) > 0
        ICpos = ICmid;
    else
        ICneg = ICmid;
    end

    if norm(ICpos-ICneg) < 1e-14*norm(ICmid)
        yShoot = ICmid;           % converged
        break
    end
end

%% 4.  Integrate the converged orbit until encounter time TE
[TE,YShoot] = ode45(@(t,y)QuadrupleRHS(t,y,G,m), ...
                    tgrid, yShoot, options);


%% 5-NEW  -- energy-error panel ---------------------------------------
% compute total energy at every stored point
E0   = totalEnergy(YShoot(1,:)',G,m);          % reference energy
Ener = arrayfun(@(k) totalEnergy(YShoot(k,:)',G,m), 1:size(YShoot,1));

relErr = abs(Ener - E0)./abs(E0);

figure(3); clf
semilogy(TE, relErr, 'k','LineWidth',1.2)
xlabel('t  [yr]'); ylabel('|\Delta E/E_0|')
title('Relative energy error along homoclinic shoot')
grid on
fprintf('max |ΔE/E| = %.2e over %.2f yr\n', max(relErr), TE(end));



%---------------------------------------------------------
% 5.  Plot outer 4-body motion up to closest approach
%---------------------------------------------------------
figure(1); clf; hold on
plot(YShoot(:,1), YShoot(:,2),'r','LineWidth',1.3)      % star A
plot(YShoot(:,4), YShoot(:,5),'r:','LineWidth',1.0)     % star B
plot(YShoot(:,7), YShoot(:,8),'b','LineWidth',1.3)      % star C
plot(YShoot(:,10),YShoot(:,11),'b:','LineWidth',1.0)    % star D
axis equal
xlabel('x [AU]'); ylabel('y [AU]');
title('Homoclinic-like outer encounter of AB and CD');

%---------------------------------------------------------
% 6.  Single closed orbit of inner binary AB for validation
%---------------------------------------------------------
optsNoE = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,orbitAB] = ode45(@(t,y)QuadrupleRHS(t,y,G,m), [0 Pab], y0, optsNoE);

figure(2); clf
plot( orbitAB(:,1)-orbitAB(:,4), orbitAB(:,2)-orbitAB(:,5), 'k-' );
axis equal
xlabel('x [AU]'); ylabel('y [AU]');
title('Inner binary AB – one period');

% ===================  nested function definitions  ===================

function dY = QuadrupleRHS(~,Y,G,m)
    dY = zeros(24,1);
    r = reshape(Y(1:12),3,4);   v = reshape(Y(13:24),3,4);
    for i = 1:4
        a = [0;0;0];
        for j = 1:4
            if i~=j
                rij = r(:,i)-r(:,j);
                a = a - G*m(j)*rij/norm(rij)^3;
            end
        end
        dY(3*(i-1)+(1:3))      = v(:,i);
        dY(12+3*(i-1)+(1:3))   = a;
    end
end

function [value,isterm,dir] = CloseApproach(~,Y,G,m)
    r  = reshape(Y(1:12),3,4);
    v  = reshape(Y(13:24),3,4);
    rAB = 0.5*(r(:,1)+r(:,2));   vAB = 0.5*(v(:,1)+v(:,2));
    rCD = 0.5*(r(:,3)+r(:,4));   vCD = 0.5*(v(:,3)+v(:,4));
    Rvec  = rAB - rCD;           Vrel = vAB - vCD;
    value = Rvec.'*Vrel / norm(Rvec);   % dR/dt
    isterm = 1; dir = +1;               % stop at first minima
end

% ---------- totalEnergy ----------------------------------------------
function E = totalEnergy(y,G,m)
% y = 24×1 state vector  [r1 … r4 v1 … v4]^T  (row or column OK)
y = y(:);
r = reshape(y(1:12) ,3,4);       % positions
v = reshape(y(13:24),3,4);       % velocities

% kinetic energy
KE = 0.5*sum( m(:)'.*sum(v.^2,1) );

% potential energy (pairwise −G m_i m_j / |r_i−r_j|)
PE = 0;
for i = 1:3
    for j = i+1:4
        rij = norm(r(:,i)-r(:,j));
        PE  = PE - G*m(i)*m(j)/rij;
    end
end
E = KE + PE;
end


end
