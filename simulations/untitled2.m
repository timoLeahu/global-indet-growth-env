function untitled2()
 IC = rand(10,2);
 hold on
 for ii = 1:length(IC(:,1))
    [~,X] = ode45(@EOM,[0 2],IC(ii,:));
    x = X(:,1);
    y = X(:,2);
    plot(x,y,'r')
 end
 xlabel('x')
 ylabel('y')
 grid
 end
 function dZ = EOM(t, z)
 dZ = zeros(2,1);
 x  = z(1);
 y  = z(2);
 dZ = [  y;...
       - x - x^2];
 end