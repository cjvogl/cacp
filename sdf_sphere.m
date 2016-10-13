function [phi, nx, ny, nz, kappa, h] = sdf_sphere(X,Y,Z)

    psi = sqrt(X.^2 + Y.^2 + Z.^2) - 1;
    cval = classify(psi);
    ind = cval > 0;
    
    phi = zeros(size(psi));
    nx = zeros(size(psi));
    ny = zeros(size(psi));
    nz = zeros(size(psi));
    kappa = zeros(size(psi));
    h = zeros(size(psi));
    
    phi(ind) = sqrt(X(ind).^2 + Y(ind).^2 + Z(ind).^2) - 1.0;
    phi(~ind) = sign(psi(~ind))*inf;
    nx(ind) = X(ind)./sqrt(X(ind).^2 + Y(ind).^2 + Z(ind).^2 + eps);
    ny(ind) = Y(ind)./sqrt(X(ind).^2 + Y(ind).^2 + Z(ind).^2 + eps);
    nz(ind) = Z(ind)./sqrt(X(ind).^2 + Y(ind).^2 + Z(ind).^2 + eps);
    kappa(ind) = 1;
    h(ind) = 1;

end