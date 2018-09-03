rad
====


A large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.

:\Esun_: incident solar spectrum

    :units: mW m-2 um-1
    :type: [2162x1 double]

:\Esky_: [2162x1 double]   incident sky spectrum (mW m-2 um-1)
:PAR: [1 double]        incident spectrally integrated PAR (moles m-2 s-1)

:fEsuno:
    [2162x1 double]   normalized spectrum of direct light (optical)
:fEskyo:
    [2162x1 double]   normalized spectrum of diffuse light (optical)
:fEsunt:
    [2162x1 double]   normalized spectrum of direct light (thermal)
:fEskyt:
    [2162x1 double]   normalized spectrum of diffuse light (thermal)
:\Eplu_:
    [61x2162 double]  upward diffuse radiation in the canopy (mW m-2 um-1)
:\Emin_:
    [61x2162 double]  downward diffuse radiation in the canopy (mW m-2 um-1)
:\Lo_:
    [2162x1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)
:\Eout_:
    [2162x1 double]   TOC upward radiation (mW m-2 um-1)
:Eouto:
    [1 double]        TOC spectrally integrated upward optical ratiation (W m-2)
:Eoutt:
    [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)

:Rnhs:
    [1 double]        net radiation (W m-2) of shaded soil
:Rnus:
    [1 double]        net radiation (W m-2) of sunlit soil
:Rnhc:
    [60x1 double]     net radiation (W m-2) of shaded leaves
:Rnuc:
    [13x36x60 double] net radiation (W m-2) of sunlit leaves
:Pnhc:
    [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
:Pnuc:
    [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
:Pnhc_Cab:
    [60x1 double]     net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
:Pnuc_Cab:
    [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
:Rnhc_PAR:
    [60x1 double]     net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
:Rnuc_PAR:
    [13x36x60 double] net PAR absorbed (W m-2) of sunlit