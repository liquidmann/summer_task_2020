x0 = 14;  %initial number of prey
y0 = 14;  %initial number of predators
z0 = 20; %initial number of food for prey

x_goal = 15; %prey goal number
B = 20; %predator restriction


fid = fopen( 'numbers.txt', 'wt' );

left = 1; %left interval border
right = 100; %right interval border
step = 0.01; %step size
mass = left:step:right;
size = length(mass);

x = zeros(1,size);
y = zeros(1,size);
z = zeros(1,size);

x(1) = x0;
y(1) = y0;
z(1) = z0;

for T1 = 0.01:0.01:1
    for T2 = 1:1:10
        for beta = 0.1:0.1:1.5
            for gamma = 0.1:0.1:1.5
                for sigma = 0.1:0.1:1.5
                    
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
                                        
                                        
                                        
                                        rest = size - size / 100;
                                        
                                        if(isnan(x(i)) || isnan(y(i)))
                                            break
                                        elseif(i > rest && abs(psi_t)/x_goal > 0.01)
                                            break
                                        elseif(abs(y(i)) > B)
                                            break
                                        elseif(i==size)
                                            fprintf( fid, '%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f\n', T1, T2, beta, gamma, sigma, x0, y0, z0, x_goal, B);

                                        end

                      end
                end
            end
        end
    end
end

fclose(fid); 