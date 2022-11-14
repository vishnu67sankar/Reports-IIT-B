function g = def_g_uf(t, eta, eta_t, eta_count_uf, tau)
    [gamma, c1, c2, M, k, xf, J, w, zeta, uf, mode_size] = initialising_variables;
    
    
    if (t>tau)           
           for iter = 1:mode_size
               uf = uf + eta_count_uf(iter)*cos(iter*pi*xf);               
           end
    end
    
    g = -2.*zeta.*w.*eta_t - (w.*w).*eta - (2*k/(gamma*M)).*(pi*J)*((abs(1/3 + uf))^0.5 - (1/3)^0.5).*sin(J*pi*xf);
    
end