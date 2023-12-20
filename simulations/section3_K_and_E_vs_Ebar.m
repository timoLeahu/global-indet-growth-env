clear; close all
% generation of two graphs of the variation of coordinates K & E
% of the fixed points P* of the model: plot K-E_bar & E-E_bar

% params
alpha = 0.1;
beta = 0.8; 
gamma = 0.58;
delta = 0.05;
epsil = 1;
eta = 1.5;
theta = 0.001;
%E_bar = 0.21;
E_temp = 0.17:0.01/2:0.6;
E_bar(1)=0.1671597;
for j=1:size(E_temp,2)
    E_bar(j+1)=E_temp(j);
end

syms x y z
for i=1:size(E_bar,2)
    eqn1 = z-(beta/(beta+epsil)) == 0;
    eqn2 = x-(alpha/(delta*theta))*y*(E_bar(i)-y) == 0;
    eqn3 = y+delta*(beta/(beta+epsil))^(beta/(1-alpha))*...
    (alpha/theta)^(alpha/(1-alpha))*y^((alpha+gamma-1)/(1-alpha))-E_bar(i) == 0;

    sol = solve([eqn1, eqn2, eqn3], [x, y, z],'ReturnConditions', true, 'Real',true);
    xSol = vpa(sol.x); % K
    ySol = vpa(sol.y); % E
    zSol = vpa(sol.z); % L
    if xSol
    K1(i)=xSol(1); % K values for fixed point P1*
    E1(i)=ySol(1); % E values for fixed point P1*
    K2(i)=xSol(2); % K values for fixed point P2*
    E2(i)=ySol(2); % E values for fixed point P2*
    end
end
figure(1) %plot K vs E_bar
hold on
p(1)=plot(E_bar(1:8),K1(1:8),'LineWidth',2,'Color','k'); % trajectory of P1*(sink)
p(2)=plot(E_bar(8:size(E_bar,2)),K1(8:size(K1,2)),'-.','LineWidth',1.5);% P1*(saddle)
p(3)=plot(E_bar(1:6),K2(1:6),':','LineWidth',1); % trajectory of P2*(saddle)
p(4)=plot(E_bar(1), K1(1)+0.01, 'Marker','o','Color','k');
p(5)=plot(0.2, 5.12116, 'Marker','*','Color','k');
xlabel('$\bar{\textbf{E}}$', 'Interpreter','latex')
ylabel('\bf K', 'Interpreter','latex')
legend(p,{'Locus of P_1^* as a sink','Locus of P_1^* as a saddle','Locus of P_2^* as a saddle',...
    'Lowest Point (LP)','Hopf bifurcation point (H)'},'location','northeast')
box on
hold off

figure(2) %plot E vs E_bar
hold on
p(1)=plot(E_bar(1:8),E1(1:8),'LineWidth',2,'Color','k');% trajectory of P1*(sink)
p(2)=plot(E_bar(8:size(E_bar,2)),E1(8:size(E1,2)),'-.','LineWidth',1.5);% P1*(saddle)
p(3)=plot(E_bar(1:6),E2(1:6),':','LineWidth',1);% trajectory of P2*(saddle)
p(4)=plot(E_bar(1), E1(1), 'Marker','o','Color','k');
p(5)=plot(0.2, 0.0137479, 'Marker','*','Color','k');
xlabel('$\bar{\textbf{E}}$', 'Interpreter','latex')
ylabel('\bf E', 'Interpreter','latex')
legend(p,{'Locus of P_1^* as a sink','Locus of P_1^* as a saddle','Locus of P_2^* as a saddle',...
    'Lowest Point (LP)','Hopf bifurcation point (H)'},'location','northeast')
box on