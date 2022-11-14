function [eta, eta_t] = RK4( count, t, delta_t, eta, eta_t )
       
       k1 = def_f(t, eta(count,:), eta_t(count,:)); 
       l1 = def_g(t, eta(count,:), eta_t(count,:)); 
       
       k2 = def_f(t + delta_t/2, eta(count,:)+ k1/2, eta_t(count,:) + l1/2); 
       l2 = def_g(t + delta_t/2, eta(count,:)+ k1/2, eta_t(count,:) + l1/2);
       
       k3 = def_f(t + delta_t/2, eta(count,:)+ k2/2, eta_t(count,:) + l2/2); 
       l3 = def_g(t + delta_t/2, eta(count,:)+ k2/2, eta_t(count,:) + l2/2);
       
       k4 = def_f(t + delta_t/2, eta(count,:)+ k3/2, eta_t(count,:) + l3/2);
       l4 = def_g(t + delta_t/2, eta(count,:)+ k3/2, eta_t(count,:) + l3/2);
       
       eta(count+1,:) = eta(count,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
       eta_t(count+1,:) = eta_t(count,:) + (1/6)*(l1 + 2*l2 + 2*l3 + l4);
       
end