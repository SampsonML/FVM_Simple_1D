function c = ANALYTIC_ADVECTDIFF(x,t,u,D,lambda,c0)
% Returns the analytic solution to an advection diffusion equation%

gamma = sqrt(1 + (4*lambda*D)/(u^2));

c = (c0 / 2) .* ( exp( (u.*(1-gamma).*x) ./ (2*D) ) .* erfc( (x - u*gamma*t ) ./ (2*sqrt(D*t) ) ) ...
                                + exp( (u.*(1+gamma).*x) ./ (2*D) ) .* erfc( (x + u*gamma*t ) ./ (2*sqrt(D*t) ) )  );

end