function Lb = Planck (wl,Tb,em)

    c1 = 1.191066e-22;
    c2 = 14388.33;
    if nargin<3
        em = ones(size(Tb));
    end
    
    Lb = em.* c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tb))-1);
    
end    