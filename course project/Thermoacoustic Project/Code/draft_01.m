clear all
close all
clc
%% Variable Definitions
% gamma = 1.4;
% c1 = 0.1;
% c2 = 0.06;
% M = 0.1;
% k = 0.5;
% xf = 0.25;
% J = (1:10);
% w = J*pi;
% zeta = (1/(2*pi))*(c1*w/w(1) + c2*(w(1)/w).^0.5);
% uf = 0.0;

[gamma, c1, c2, M, k, xf, J, w, zeta, uf] = initialising_variables;
Tmax = 40; %10 seconds
window = 100001;
t_window = linspace(0,Tmax,window);
delta_t = (t_window(end) - t_window(1))/(window-1);
tau = 0.45*pi; 

%% Variable Definitions
% for i = 1:10
%    eta_0(1) = rand(); 
%    eta_t1(i) = rand(); 
% end

% eta_0(1) = 3.0; // testing RK4
% eta_t1(1) = 1.0; // testing RK4
eta_0(1) = 0.15; %0.15;
eta_0(2:10) = 0.0;
eta_t1(1) = 0.0;
eta_t1(2:10) = 0.0;

eta = zeros(window,10); 
eta_t = zeros(window,10);
eta(1,:) = eta_0;
eta_t(1,:) = eta_t1;

%% Range Kutta
t = 0;
count_uf = 0;
for count = 1:window-1
       
%        [eta, eta_t] = RK4(count, t, delta_t, eta(count,:), eta_t(count,:)); 
       %count_f = max(ceil((t - tau)/delta_t),1);
       
       k1 = delta_t*def_f(t, eta(count,:), eta_t(count,:));        
       l1 = delta_t*def_g(t, eta(count,:), eta_t(count,:)); 
       
       k2 = delta_t*def_f(t + delta_t/2, eta(count,:) + k1/2, eta_t(count,:) + l1/2); 
       l2 = delta_t*def_g(t + delta_t/2, eta(count,:) + k1/2, eta_t(count,:) + l1/2);
       
       k3 = delta_t*def_f(t + delta_t/2, eta(count,:) + k2/2, eta_t(count,:) + l2/2); 
       l3 = delta_t*def_g(t + delta_t/2, eta(count,:) + k2/2, eta_t(count,:) + l2/2);
       
       k4 = delta_t*def_f(t + delta_t/2, eta(count,:) + k3/2, eta_t(count,:) + l3/2);
       l4 = delta_t*def_g(t + delta_t/2, eta(count,:) + k3/2, eta_t(count,:) + l3/2);
       
       eta(count+1,:) = eta(count,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
       eta_t(count+1,:) = eta_t(count,:) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
    
    
    t = t + delta_t;

end

%plot(t_window(1:5), eta(1:5,1));
%plot(t_window(:), real(eta(:,1)));


u = zeros(window,1);
freq = sin(J*pi*xf);

for iter = 1:window
   u(iter,1) = u(iter,1) +  sum(eta(iter,:).*freq);
end

figure()
plot(t_window(:), real(u(:,1)));
xlabel('t (s)');
ylabel('u_{hat} (m/s)');
title('Velocity Fluctuations vs. Time');