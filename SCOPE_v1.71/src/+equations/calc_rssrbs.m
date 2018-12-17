function [rss,rbs] = calc_rssrbs(SMC,LAI,rbs)

rss            = 11.2*exp(42*(0.22-SMC));
rbs            = rbs*LAI/3.3;