function [GAM]  = Soil_Inertia0(cs,rhos,lambdas)
% soil thermal intertia 
GAM             = sqrt(cs*rhos*lambdas);                            % soil thermal intertia