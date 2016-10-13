function [phi, nx, ny, kappa] = sdf_clover(X,Y)

    global alpha

    theta = atan2(Y,X);
    alpha = 0.25;
    
    psi = X.^2 + Y.^2 - (1 + alpha*cos(4*theta-pi)).^2;
    cval = classify(psi);
    ind = find(cval > 0);

    phi = sign(psi)*inf;
    nx = zeros(size(psi));
    ny = zeros(size(psi));
    kappa = zeros(size(psi));
    
    for iter = 1:length(ind)
        k = ind(iter);
        objf = @(val)(sqrt((X(k) - (1.0 + alpha*cos(4*val-pi)).*cos(val)).^2 ...
                          +(Y(k) - (1.0 + alpha*cos(4*val-pi)).*sin(val)).^2));
        options = optimset('TolFun',1e-15);
        [theta_star,aphi2,eflag] = fminsearch(objf,theta(k),options);
        if (eflag ~= 1)
            error('Distance minimization failed')
        end
        f = 1.0 + alpha*cos(4*theta_star-pi);
        fp = -4.0*alpha*sin(4*theta_star-pi);
        fpp = -16.0*alpha*cos(4*theta_star-pi);
        phi(k) = aphi2;
        nx(k) = -fp*(-sin(theta_star)) + f*cos(theta_star);
        nx(k) = nx(k)/sqrt(fp^2 + f^2);
        ny(k) = -fp*cos(theta_star) + f*sin(theta_star);
        ny(k) = ny(k)/sqrt(fp^2 + f^2);
        kappa(k) = (f^2 - fpp*f + 2*fp^2)/(fp^2 + f^2)^1.5;
    end

    phi(ind) = phi(ind).*sign(psi(ind));
end    
    
    
