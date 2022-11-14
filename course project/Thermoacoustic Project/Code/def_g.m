function g = def_g(t, eta, eta_t)
    [gamma, c1, c2, M, k, xf, J, w, zeta, uf] = initialising_variables;
    tau = 0.45*pi;
    
    if (t>tau)
           
           for iter = 1:10
               uf = uf + eta(iter)*cos(iter*pi*xf);               
           end
    end
    
    g = -2.*zeta.*w.*eta_t - (w.*w).*eta - (2*k/(gamma*M)).*(pi*J)*((1/3 + uf)^0.5 - (1/3)^0.5).*sin(J*pi*xf);
    
end