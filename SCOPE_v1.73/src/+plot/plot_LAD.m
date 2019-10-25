function plot_LAD(a, b, name)
    prob = equations.leafangles(a, b);
    % angles = [10:10:80, 82:2:88 90];
    % angles_mirrowed = [92:2:98, 100:10:180];
    prob_agg = [prob(1:8); sum(prob(9:end))];
    polarhistogram('BinEdges', deg2rad([0 10:10:180]), ...
        'BinCounts',[prob_agg; flip(prob_agg)], 'FaceColor', 'g')
    thetalim([0, 180])
    set(gca,'RTickLabel',[])
    title(name)
    legend(sprintf('LIDFa = % .2f\nLIDFb = % .2f', a, b))
end