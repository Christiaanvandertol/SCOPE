function Lb = Planck(wl,Tb,em)

    c1 = 1.191066e-22;  % 2 * h * c ^ 2 * (1e-6)  to make the resulting radiance per um
    c2 = 14388.33;      % h * c / k  * (1e6)  to make it um * K
    if nargin<3
        em = ones(size(Tb));
    end
    
    Lb = em.* c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tb))-1);
    
end    