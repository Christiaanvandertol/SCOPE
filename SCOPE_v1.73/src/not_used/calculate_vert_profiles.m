function profiles = calculate_vert_profiles(profiles, canopy)
% this function is incomplete and apparently never called

profiles.etah = Fh;
profiles.etau = Fu;

[Hcu1d  ]           = equations.meanleaf(canopy,Hcu,          'angles');   % [nli,nlo,nl]      mean sens heat sunlit leaves
[lEcu1d ]           = equations.meanleaf(canopy,lEcu,         'angles');   % [nli,nlo,nl]      mean latent sunlit leaves
[Au1d   ]           = equations.meanleaf(canopy,Au,           'angles');   % [nli,nlo,nl]      mean phots sunlit leaves
[Fu_Pn1d]           = equations.meanleaf(canopy,Fu.*Pinu_Cab, 'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
[qEuL   ]           = equations.meanleaf(canopy,qEu,          'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
[Pnu1d  ]           = equations.meanleaf(canopy,Pinu,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
[Pnu1d_Cab  ]       = equations.meanleaf(canopy,Pinu_Cab,     'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
[Rnu1d  ]           = equations.meanleaf(canopy,Rncu,         'angles');   % [nli,nlo,nl]      mean net PAR sunlit leaves
[Tcu1d  ]           = equations.meanleaf(canopy,Tcu,          'angles');   % [nli,nlo,nl]      mean temp sunlit leaves

profiles.Tchave     = mean(Tch);                                           % [1]     mean temp shaded leaves
profiles.Tch        = Tch;                                                 % [nl]
profiles.Tcu1d      = Tcu1d;                                               % [nl]
profiles.Tc1d       = (1-Ps(1:nl)).*Tch       + Ps(1:nl).*(Tcu1d);         % [nl]    mean temp leaves, per layer
profiles.Hc1d       = (1-Ps(1:nl)).*Hch       + Ps(1:nl).*(Hcu1d);         % [nl]    mean sens heat leaves, per layer
profiles.lEc1d      = (1-Ps(1:nl)).*lEch      + Ps(1:nl).*(lEcu1d);        % [nl]    mean latent heat leaves, per layer
profiles.A1d        = (1-Ps(1:nl)).*Ah        + Ps(1:nl).*(Au1d);          % [nl]    mean photos leaves, per layer
profiles.F_Pn1d     = ((1-Ps(1:nl)).*Fh.*Pinh_Cab + Ps(1:nl).*(Fu_Pn1d));  % [nl]    mean fluor leaves, per layer
profiles.qE         = ((1-Ps(1:nl)).*qEh      + Ps(1:nl).*(qEuL));         % [nl]    mean fluor leaves, per layer
profiles.Pn1d       = ((1-Ps(1:nl)).*Pinh     + Ps(1:nl).*(Pnu1d));        % [nl]    mean aPAR leaves, per layer
profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pinh_Cab + Ps(1:nl).*(Pnu1d_Cab));    % [nl]    mean aPAR_byCab leaves, per layer
profiles.Rn1d       = ((1-Ps(1:nl)).*Rnch     + Ps(1:nl).*(Rnu1d));        % [nl]    mean net radiation leaves, per layer
end