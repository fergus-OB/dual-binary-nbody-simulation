function [G,m,y0,Pab] = hd74438_setup
% HD74438 orbital and mass constants in code units
% Returns:
%   G   – gravitational constant (AU^3 M⊙⁻¹ yr⁻²)
%   m   – 1×4 vector of stellar masses (M⊙)
%   y0  – 24×1 initial-state vector [r1 … r4 v1 … v4]^T
%   Pab – inner-binary AB period (code units)

% -------- constants --------
G  = 39.47841760435743;
m  = [1.65 1.65 1.65 1.65];     % A, B, C, D
Pab = 4.35/365.25;  Pcd = 20.5/365.25;  Pout = 5.7;
eab = 0.15;         ecd = 0.36;         eout = 0.46;

% -------- semi-major axes (AU) --------
aab  = (G*(m(1)+m(2))*Pab^2 /(4*pi^2))^(1/3);
acd  = (G*(m(3)+m(4))*Pcd^2 /(4*pi^2))^(1/3);
aout = (G*sum(m)*Pout^2/(4*pi^2))^(1/3);

% -------- Kepler → Cartesian (true anomaly f=0) --------
rAB  = aab *(1-eab);     vAB  = sqrt(G*(m(1)+m(2))*(1+eab)/aab /(1-eab));
rCD  = acd *(1-ecd);     vCD  = sqrt(G*(m(3)+m(4))*(1+ecd)/acd /(1-ecd));
rOUT = aout*(1-eout);    vOUT = sqrt(G*sum(m)*(1+eout)/aout/(1-eout));

% positions (planar)
r1 = [+m(2)/(m(1)+m(2))*rAB, 0, 0];
r2 = [-m(1)/(m(1)+m(2))*rAB, 0, 0];
r3 = [+rOUT - m(4)/(m(3)+m(4))*rCD, 0, 0];
r4 = [+rOUT + m(3)/(m(3)+m(4))*rCD, 0, 0];

% velocities
v1 = [0,+m(2)/(m(1)+m(2))*vAB, 0];
v2 = [0,-m(1)/(m(1)+m(2))*vAB, 0];
v3 = [0,-vOUT - m(4)/(m(3)+m(4))*vCD, 0];
v4 = [0,-vOUT + m(3)/(m(3)+m(4))*vCD, 0];

% shift to barycentre
Rcm = (m(1)*r1 + m(2)*r2 + m(3)*r3 + m(4)*r4)/sum(m);
Vcm = (m(1)*v1 + m(2)*v2 + m(3)*v3 + m(4)*v4)/sum(m);

r1=r1-Rcm; r2=r2-Rcm; r3=r3-Rcm; r4=r4-Rcm;
v1=v1-Vcm; v2=v2-Vcm; v3=v3-Vcm; v4=v4-Vcm;

% 24-element state vector
y0 = [r1 r2 r3 r4 v1 v2 v3 v4]';
end
