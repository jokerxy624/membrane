% Numerical solution of the capillary attraction-driven kinematic equation
% Definition and description ([]: physical units, @: experimental value)

clear
clc

% Variables
x = ;       % Initial distance between the two floating objects, [m]
t = ;       % Duration of observation, [s]
dt = 0.0001;       % Temporal resolution of iteration interval, [s]
r1 = ;       % Radius of spherical floating object 1, [m]
K = ;       % Linear diversity factor
r2 = K*r1;       % Radius of spherical floating object 2, [m]
p = ;       % @Density of the material, [kg m-3]
m1 = 4.19*(r1^3)*p;
m2 = 4.19*(r2^3)*p;       % Masses of floating object 1 and 2, [kg]

% Constants
f = 0.0167676;       % Coefficient of viscous force
G = 1.45*(p^4)*(r1^3)*(r2^3);       % Equivalent attraction coefficient
B = 0.05;       % Correction factor before G, the portion of effective attraction
a = 2.68E-3;       % Capillary length, [m]

% Initial settings
X = [];
v1 = 0;
v2 = 0;       % Both of the floating objects are still, namely initial velocities are 0
a1 = 0;
a2 = 0;       % Initial accelerated velocities of both of the floating objects are 0
x1 = 0;
x2 = 0;       % Initial displacement of both of the floating objects are 0
V1 = [];
V2 = [];

for i = 0:t/dt
     
     v1 = v1+a1*dt;
     v2 = v2+a2*dt;
     f1 = f*r1*v1;       % Viscous resistance, [N]
     f2 = f*r2*v2;
     A = B*G*besselk(1,x/a);       % Attractive force, [N]
     a1 = (A-f1)/m1;
     a2 = (A-f2)/m2;
     x = x-v1*dt-v2*dt;       % Variable distance, [m]
     X = [X,x];
     V1 = [V1,v1];
     V2 = [V2,-v2];
     x1 = x1+v1*dt;
     x2 = x2+v2*dt;
     if r1+r2 > x
         break       % The distance can¡¯t be less than (r1+r2) in fact
     end
end

T = 0:dt:i*dt;
plot(T,real(X));
plot(T,V1);
hold on;
plot(T,V2);
