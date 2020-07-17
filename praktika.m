function praktika(T1, T2, beta, gamma, sigma, x0, y0, z0, x_goal, B)  







left = 1;%left interval border
right = 100;%right interval border
step = 0.01;%step size
mass = left:step:right;
size = length(mass);
x = zeros(1,size);
y = zeros(1,size);
z = zeros(1,size);



x(1) = x0;   %initial number of prey
y(1) = y0;  %initial number of predators
z(1) = z0;  %initial number of food for prey






fid = fopen( 'results.txt', 'wt' );

for i = 2:size
    
    f1 = z(i-1) * x(i-1) - beta  * x(i-1) * y(i-1);
    f2 = -gamma * y(i-1) + sigma * x(i-1) * y(i-1);
    df2 = - gamma * f2 + sigma * (f1 * y(i-1) + x(i-1) * f2);
    
    psi_t = x(i-1) - x_goal;
    
    P = 4 * ((exp(psi_t) + exp(-psi_t))^(-2));
    dP = -2 * f1 * P * (exp(2 * psi_t) + exp(-2 .* psi_t));
   
    psi_0 = y(i-1) + B * tanh(psi_t);
    dpsi_0 = f2 + P * f1;
    
    dfi = (-(T2^(-1) * dpsi_0 + df2) * B * P * x(i-1) + B * (dP * x(i-1) + P * f1) * (T2^(-1) * psi_0 + f2)) / ((B * P * x(i-1))^2) + beta * f2;
    
    fi = - T2^(-1) + x_goal / (T2 * x(i-1)) + beta .* y(i-1);
    psi_1 = z(i-1) - fi;
    U = - T1^(-1) * psi_1 + dfi;
    
    
    x(i) = x(i-1) + step * f1;
    y(i) = y(i-1) + step * f2;
    z(i) = z(i-1) + step * U;
    
    
end

fclose(fid); 


mass_plot = 1:size;
hold on;
plot(mass_plot, x, 'b', mass_plot, y, 'r', mass_plot, z, 'g');
legend('Prey', 'Predator', 'Food');
title('Predator-prey with food');
ylabel('X,Y,Z');
grid on;
hold off;
end
