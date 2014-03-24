function myBoxplot(X, colors, markers)
% Written by Emtiyaz,
% May 31, 2010

  [nInstances, nVars] = size(X);

  med = median(X,1);

  for m = 1:nVars
    hold on
    errQ1 = quantile(X(:,m),0.25);
    errQ2 = quantile(X(:,m),0.75);
    line([m, m], [errQ1 errQ2], 'linewidth',3,'color',colors{m});
    plot(m, med(m), 'color', colors{m}, 'marker', markers{m}, 'linewidth',3, 'markersize',12, 'markerfacecolor',[1 1 1]);
  end
  xlim([0.5 nVars+0.5])

