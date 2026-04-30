
% -------------- hd74438_setup.m --------------------
% Constants (AU, Msun, yr)
G  = 39.47841760435743;        % AU^3 Msun^-1 yr^-2% AU^3 Msun^-1 yr^-2  m  = [1.65 1.65 1.65 1.65];   % old
% ↓ new (Merle et al. 2022, Nat. Astron. 6, 681)
m  = [1.70 1.54 0.96 0.87];
Pab = 4.35/365.25;  Pcd = 20.5/365.25; Pout = 5.7;
eab = 0.15;         ecd = 0.36;         eout = 0.46;

% Semi‑major axes
aab  = (G*(m(1)+m(2))*Pab^2/(4*pi^2))^(1/3);
acd  = (G*(m(3)+m(4))*Pcd^2/(4*pi^2))^(1/3);
aout = (G*sum(m)*Pout^2/(4*pi^2))^(1/3);

% Keplerian-to-Cartesian for each binary (true anomaly f=0)
rAB  = aab*(1-eab);   vAB  = sqrt(G*(m(1)+m(2))*(1+eab)/aab/(1-eab));
rCD  = acd*(1-ecd);   vCD  = sqrt(G*(m(3)+m(4))*(1+ecd)/acd/(1-ecd));
rOUT = aout*(1-eout); vOUT = sqrt(G*sum(m)*(1+eout)/aout/(1-eout));

% Position vectors (planar: AB along +x, CD along -x, outer along +x)
r1 = [+m(2)/(m(1)+m(2))*rAB, 0, 0];
r2 = [-m(1)/(m(1)+m(2))*rAB, 0, 0];
r3 = [+rOUT - m(4)/(m(3)+m(4))*rCD, 0, 0];
r4 = [+rOUT + m(3)/(m(3)+m(4))*rCD, 0, 0];

% Velocities (y‑direction for inner binaries, +y for AB, -y for CD)
v1 = [0, +m(2)/(m(1)+m(2))*vAB, 0];
v2 = [0, -m(1)/(m(1)+m(2))*vAB, 0];
v3 = [0, -vOUT - m(4)/(m(3)+m(4))*vCD, 0];
v4 = [0, -vOUT + m(3)/(m(3)+m(4))*vCD, 0];
% ----------------------------------------------------------
% shift to barycentric coordinates (total P = 0 , R = 0)
Rcm = (m(1)*r1 + m(2)*r2 + m(3)*r3 + m(4)*r4)/sum(m);
Vcm = (m(1)*v1 + m(2)*v2 + m(3)*v3 + m(4)*v4)/sum(m);

r1 = r1 - Rcm;  r2 = r2 - Rcm;  r3 = r3 - Rcm;  r4 = r4 - Rcm;
v1 = v1 - Vcm;  v2 = v2 - Vcm;  v3 = v3 - Vcm;  v4 = v4 - Vcm;
% ----------------------------------------------------------

y0 = [r1 r2 r3 r4 v1 v2 v3 v4]';  % 24‑vector
% --------- choose the scenario -------------
outerFull   = true;   % 6 yr
singleInner = false;   % ~one inner orbit
longRun     = false;    % 200 yr
% -------------------------------------------

if outerFull
    tspan = [0 6];
elseif singleInner
    tspan = [0 0.012];          % 4.35 d  ≃ 0.012 yr   % change manually for CD
elseif longRun
    tspan = [0 200];
end

                        % 200 yr ≈ 35 outer periods
opts  = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,Y] = ode45(@rhs_4body, tspan, y0, opts, G, m);
% ---- Outer barycentre trajectories ---------------------------------
rAB_cm = (m(2)*Y(:,1:3) + m(1)*Y(:,4:6))/(m(1)+m(2));
rCD_cm = (m(4)*Y(:,7:9) + m(3)*Y(:,10:12))/(m(3)+m(4));

figure; hold on; axis equal
plot(rAB_cm(:,1), rAB_cm(:,2), 'r');
plot(rCD_cm(:,1), rCD_cm(:,2), 'b');
plot(0,0,'k+','MarkerSize',6);
xlabel('x [AU]'); ylabel('y [AU]');
title('Outer orbit of HD 74438 (barycentric)');
legend('Binary AB CM','Binary CD CM','System CM', 'Location','Best');

% ----------- Zoom: Inner binary AB (relative coordinates) -------------
r1rel = Y(:,1:3)  - rAB_cm;     % star A trajectory minus AB CM
r2rel = Y(:,4:6)  - rAB_cm;     % star B trajectory minus AB CM
figure; axis equal
plot(r1rel(:,1), r1rel(:,2), 'r'); hold on
plot(r2rel(:,1), r2rel(:,2), 'r');
xlim(0.08*[-1 1]); ylim(0.08*[-1 1]);
xlabel('x [AU]'); ylabel('y [AU]');
title('Inner binary AB (a ≈ 0.066 AU)');

% ----------- Zoom: Inner binary CD (relative coordinates) -------------
r3rel = Y(:,7:9)  - rCD_cm;     % star C relative to CD CM
r4rel = Y(:,10:12)- rCD_cm;     % star D relative to CD CM
figure; axis equal
plot(r3rel(:,1), r3rel(:,2), 'b'); hold on
plot(r4rel(:,1), r4rel(:,2), 'b');
xlim(0.25*[-1 1]); ylim(0.25*[-1 1]);
xlabel('x [AU]'); ylabel('y [AU]');
title('Inner binary CD (a ≈ 0.17 AU)');






% -------------- rhs_4body.m -------------------------
function dy = rhs_4body(~,y,G,m)
dy = zeros(24,1);
r = reshape(y(1:12),3,4);
v = reshape(y(13:24),3,4);
for i=1:4
   acc = [0;0;0];
   for j=1:4
       if i~=j
           rij = r(:,i)-r(:,j);
           acc = acc - G*m(j)*rij/norm(rij)^3;
       end
   end
   dy(3*(i-1)+1:3*i)     = v(:,i);
   dy(12+3*(i-1)+1:12+3*i)= acc;
end
end
