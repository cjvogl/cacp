function [phi, nx, ny, kappa, theta_cp] = sdf_circle(X,Y)

    psi = X.^2 + Y.^2 - 1;
    cval = classify(psi);
    ind = cval > 0;
    
    phi = zeros(size(psi));
    nx = zeros(size(psi));
    ny = zeros(size(psi));
    kappa = zeros(size(psi));
    theta_cp = zeros(size(psi));
    
    phi(ind) = sqrt(X(ind).^2 + Y(ind).^2) - 1.0;
    phi(~ind) = sign(psi(~ind))*inf;
    nx(ind) = X(ind)./sqrt(X(ind).^2 + Y(ind).^2 + eps);
    ny(ind) = Y(ind)./sqrt(X(ind).^2 + Y(ind).^2 + eps);
    kappa(ind) = 1;
    theta_cp(ind) = atan2(Y(ind),X(ind));
    
end