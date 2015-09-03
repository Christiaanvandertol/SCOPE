function ftts = calc_ftts(P,ts)
ts = sin(ts/180*pi);
ftts = P(1)*ts.^4 + P(2)*ts.^3 + P(3)*ts.^2 + P(4)*ts + P(5);
