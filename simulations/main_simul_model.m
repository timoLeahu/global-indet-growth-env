clear; clc; close all;

%% Specify parameter set
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
E_bar = 0.21;

N=200;

%% specify the region of the plot for vector plot
% x==K, y==E, z==L corresp in our model
[K, E, L] = meshgrid(1:1:14, 0.01:0.01:0.06, 0.2:0.1:0.8);
f = @(x) epsil*((1-x).^((epsil-eta*(1+epsil))/eta)).*(x.^(1-beta))/beta; 
diff_f = @(x) -(epsil*((1-x).^((epsil-eta*(1+epsil))/eta))...
    .*(((beta+epsil)*eta-epsil).*x+(1-beta)*eta))./(beta*eta*(x-1).*(x.^beta));

Kdot = (K.^alpha).*(L.^(beta-1)).*(E.^gamma).*((beta+epsil)*L-beta)/epsil; %Note the use of .* and .^
Edot = E.*(E_bar-E)-delta*(K.^alpha).*(L.^beta).*(E.^gamma);
Ldot = (f(L)./diff_f(L)).*(alpha*(K.^(alpha-1)).*(L.^(beta-1)).*(E.^gamma)...
    .*((beta+epsil)*L-beta)/epsil + gamma*(E_bar-E-delta*(K.^alpha).*(L.^beta).*(E.^(gamma-1)))...
    +(theta-alpha*(K.^(alpha-1)).*(L.^beta).*(E.^gamma))/eta);
%first plot the vector plot with quiver
figure
quiver3(K,E,L,Kdot,Edot,Ldot)

%% system of equations Y=[ Y(1)=K, Y(2)=E, Y(3)=L ]
g = @(t,Y) [(Y(1)^alpha)*(Y(3)^(beta-1))*(Y(2)^gamma)*((beta+epsil)*Y(3)-beta)/epsil;
    Y(2).*(E_bar-Y(2))-delta*(Y(1).^alpha).*(Y(3).^beta).*(Y(2).^gamma);
    (f(Y(3))./diff_f(Y(3))).*(alpha*(Y(1).^(alpha-1)).*(Y(3).^(beta-1)).*(Y(2).^gamma)...
    .*((beta+epsil)*Y(3)-beta)/epsil + gamma*(E_bar-Y(2)-delta*(Y(1).^alpha).*(Y(3).^beta).*(Y(2).^(gamma-1)))...
    +(theta-alpha*(Y(1).^(alpha-1)).*(Y(3).^beta).*(Y(2).^gamma))/eta)]; 
hold on
%calculate the phase trajectories for different initial conditions
%for y0=5:10:50
[ts, ys] = ode45(g,[0 5000], [8 0.035 0.34]);
% plot of closed loop phase trajectories
plot(ys(:,1), ys(:,2))
%end
% x1 = linspace(0,20);
% y1 = x1/0.106;
% plot(x1,y1,'b')
% plot(0,15,'*')
% plot(740,50,'*')
% hold off
xlabel('K')
ylabel('E')
zlabel('L')