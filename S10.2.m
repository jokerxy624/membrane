% Collision energy of two Ag Brownian clusters driven by capillary attraction
% Definition and description ([]: physical units, @: experimental value)

clear
clc

% Variables
x0 = ;       % Initial distance between two Brownian clusters, [m]
t = ;        % Observation time, [s]
dt = ;        % Temporal resolution of iteration interval, [s]
r = ;      % @Radius of the silver particles as the building unit, [m]
K = ;       % Maximum voluminal diversity factor
E = zeros(K,K);        % Collision energy, total kinetic energies

% Constants
p = 10490;      % Density of metallic silver, [kg m-3]
M = 4.19*(r^3)*p;      % Mass of the particles, [kg]
f = 0.0167676;       % Coefficient of viscous force
B = 0.05;        % Correction factor before G, the portion of effective attraction
a = 2.68E-3;       % Capillary length, [m]

for n = 1:K
    for m = n:K
        x = x0;        % Stated initial distance, [m]
        v1 = 0;
        v2 = 0;         % Initial velocities of both silver Brownian clusters are 0 
        a1 = 0;
        a2 = 0;        % Initial accelerated velocities of both silver Brownian clusters are 0
        r1 = n^(1/2)*r;
        r2 = m^(1/2)*r;        % Equivalent radii, [m]
        G = 1.45*(p^4)*(r1^3)*(r2^3);       % Equivalent attraction coefficient
        m1 = n*M; 
        m2 = m*M;            % Masses of the two Brownian clusters, [kg]
        for i = 0:t/dt
           v1 = v1+a1*dt;
           v2 = v2+a2*dt;
           f1 = f*r1*v1;         % Viscous resistance, [N]
           f2 = f*r2*v2;
           A = B*G*besselk(1,x/a);       % Attractive force, [N]
           a1 = (A-f1)/m1;
           a2 = (A-f2)/m2;
           x = x-v1*dt-v2*dt;       % Variable distance, [m]
           if r1+r2 > x 
              E(m,n) = n*0.5*M*v1^2+m*0.5*M*v2^2;    % Calculate the collision energy once the two Brownian clusters are detected to touch each other
              E(n,m) = E(m,n);
              break
           end
        end  
    end
end    
    
x=1:1:K;
y=1:1:K;
surf(x,y,E);