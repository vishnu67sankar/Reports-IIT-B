clear all
close all
clc
%% Variable Initialisation

[gamma, c1, c2, M, k, xf, J, w, zeta, uf, mode_size] = initialising_variables;
Tmax = 200; 
window = 100001; % After experimentation
t_window = linspace(0,Tmax,window); % window time
delta_t = (t_window(end) - t_window(1))/(window-1); % h = delta_t
tau = 2*pi; % value of Tau from published paper

%% Variable Definitions (To use only for generating IC)
% for i = 1:mode_size
%    eta_0(1) = rand(); 
%    eta_t1(i) = rand(); 
% end

%% To be used for validating with the paper
eta_0(1) = 0.15;  % testing RK4 on published result
eta_t1(1) = 0.0; % testing RK4 on published result
%eta_t is the time derivative of eta

eta = zeros(window,mode_size); % stores eta at all times
eta_t = zeros(window,mode_size); % stores time derivative of eta at all times
eta(1,:) = eta_0;
eta_t(1,:) = eta_t1;

%% Range Kutta fourth order
t = 0;                  % time varying dynamically
count_uf = 0;
for count = 1:window-1
    
       count_uf = max(ceil((t - tau)/delta_t + eps),1);
       
       k1 = delta_t*def_f(t, eta(count,:), eta_t(count,:));        
       l1 = delta_t*def_g_uf(t, eta(count,:), eta_t(count,:), eta(count_uf,:), tau); 
       
       k2 = delta_t*def_f(t + delta_t/2, eta(count,:) + k1/2, eta_t(count,:) + l1/2); 
       l2 = delta_t*def_g_uf(t + delta_t/2, eta(count,:) + k1/2, eta_t(count,:) + l1/2, eta(count_uf,:), tau);
       
       k3 = delta_t*def_f(t + delta_t/2, eta(count,:) + k2, eta_t(count,:) + l2); 
       l3 = delta_t*def_g_uf(t + delta_t/2, eta(count,:) + k2, eta_t(count,:) + l2, eta(count_uf,:), tau);
       
       k4 = delta_t*def_f(t + delta_t, eta(count,:) + k3, eta_t(count,:) + l3);
       l4 = delta_t*def_g_uf(t + delta_t, eta(count,:) + k3, eta_t(count,:) + l3, eta(count_uf,:), tau);
       
       eta(count+1,:) = eta(count,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
       eta_t(count+1,:) = eta_t(count,:) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
        
       t = t + delta_t;

end

%% Evaluating u and P
u = zeros(window,1);
p = zeros(window,1);
sinf = sin(J*pi*xf);
cosf = sin(J*pi*xf);

for iter = 1:window
   u(iter,1) = u(iter,1) +  sum(eta(iter,:).*cosf);
   p(iter,1) = p(iter,1) +  sum(eta_t(iter,:).*sinf);
end

%% Plotting velocity perturbation
figure()
plot(t_window(:), real(u(:,1)));
xlabel('t (s)');
ylabel('u_{hat}');
title('Velocity Fluctuations vs. Time');

%% Plotting pressure perturbation
figure()
plot(t_window(:), real(p(:,1)));
xlabel('t (s)');
ylabel('p_{hat} ');
title('Pressure Fluctuations vs. Time');
