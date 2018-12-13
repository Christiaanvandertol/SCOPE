Brief history of the model
===========================


The SCOPE model has been developed between 2006 and 2009 by Wout Verhof, Joris Timmermans, Christiaan van der Tol, Anne Verhoef and Bob Su. The idea of the model was to develop a simulator for  hyperspectral VNIR observations, the surface energy budget and photosynthesis. Chlorophyll fluorescence has been part of the model. It was originally the idea to develop a 3-D radiative transfer scheme, but this idea was (temporally) abandoned, and 1-D remained a 1-D vertical model. This had the advantage that the well-known SAIL model could be used as a basis, which is easily invertible, does not require many parameters, is computationally efficient and sufficient in many cases.

The key elements of the model have been the extension to the thermal domain (Joris Timmermans) and the radiative transfer of fluorescence (Wout Verhoef), the simulation of sensible, latent and ground heat flux, stomatal opening and photosynthesis (Christiaan van der Tol) and an aerodynamic resistance scheme (Anne Verhoef). Model inversion tools are not also available, see for example Van der Tol et al., (2016 :cite:`VanderTol2016`). There have been several updates since the published version of the model (version 1.21) in 2009. Other people have contributed to the model development as well, including Ari Kornfeld, Joe Berry, Federico Magnani (mainly the biochemical part, but also other parts), and many users provided useful feedback and suggestions (see :ref:`acknowledgements:acknowledgements`).

:Model description: Van der Tol et al., 2009 :cite:`VanderTol2009`

:Biochemical routine: Van der Tol et al., 2014 :cite:`vandertol2014`

:Leaf radiative transfer scheme: Vilfan et al., 2016 :cite:`Vilfan2016`

:Model inversion: Van der Tol et al., 2016 :cite:`VanderTol2016`
