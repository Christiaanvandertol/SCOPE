function Rli = estimate_rli(Ta, ea, constants)
     % Christiaan van der Tol and Gabriel Norberto Parodi (January 18th 2012). 
     % Guidelines for Remote Sensing of Evapotranspiration, Evapotranspiration - 
     % Remote Sensing and Modeling, Ayse Irmak, IntechOpen, DOI: 10.5772/18582. 
     % Available from: https://www.intechopen.com/books/evapotranspiration-remote-sensing-and-modeling/guidelines-for-remote-sensing-of-evapotranspiration
     if Ta < 200
         Ta = Ta + constants.C2K;
     end
     emissivity = 1.24 .* ((ea ./ Ta) .^ (1/7));
     Rli = constants.sigmaSB .* emissivity .* (Ta .^ 4);
end
