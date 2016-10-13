function [X,Y,Ucp,Ucacp,Uex] = shiftp_equation(M, type) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function to solve an elliptic PDE in R^2 using the Closest Point and 
%   Curvature-Augmented Closest Point methods.
%
%   M       - number of grid cells in one direction
%   type    - specify surface type: 1 for a circle, 2 for a clover
%
%   (X,Y)   - meshgrid for plotting purposes
%   Ucp     - Closest Point method solution
%   Ucacp   - Curvature-Augmented Closest Point solution
%   Uex     - Exact solution (or highly accurate approximatation if clover)
%
%   See README.md for link to paper with more information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if (type == 1)
        f = @(x,y)(RHS_circle(x,y));
        sdf = @(x,y)(sdf_circle(x,y));
    elseif (type == 2)
        % need to specify f here
        f = @(x,y)(RHS_clover(x,y));
        sdf = @(x,y)(sdf_clover(x,y));
    elseif (nargin ~= 4)
        error('Need to specify type 1 for circle or 2 for clover')
    end
    
    %initialize grid
    [X,Y] = meshgrid(linspace(-2,2,M+1),linspace(-2,2,M+1));
    dx = X(1,2)-X(1,1);
    fprintf('\ndx = %f\n',dx);

    [phi, nx, ny, kappa] = sdf(X,Y);
     
    % Classify interface grid points
    cval = classify(phi);
    
    % Build system
    ind = find(cval > 0);
    L = length(ind);
    
    % Create mapping from (i,j) to linear system index
    indmapping = zeros(size(X));
    for iter = 1:L
        [i,j] = ind2sub(size(X),ind(iter));
        indmapping(i,j) = iter;
    end
    
    Eh1 = zeros(L);
    Eh3 = zeros(L);
    Lh = zeros(L);
    Acacp = zeros(L);
    bcp = zeros(L,1);
    bcacp = zeros(L,1);
    for iter=1:L
        k = ind(iter);
        [i,j] = ind2sub(size(X),k);

        % Interpolation matrices
        xcp = X(k) - phi(k)*nx(k);
        ycp = Y(k) - phi(k)*ny(k);
        icp = find(Y(:,1) > ycp,1) - 1;
        jcp = find(X(1,:) > xcp,1) - 1;
        if (abs(xcp - X(icp,jcp)) < 1e-10 && abs(ycp - Y(icp,jcp)) < 1e-10)
            l = indmapping(icp,jcp);
            Eh1(iter,l) = 1.0;
            Eh3(iter,l) = 1.0;
        elseif (abs(xcp - X(icp,jcp+1)) < 1e-10 && abs(ycp - Y(icp,jcp+1)) < 1e-10)
            l = indmapping(icp,jcp+1);
            Eh1(iter,l) = 1.0;
            Eh3(iter,l) = 1.0;
        elseif (abs(xcp - X(icp+1,jcp)) < 1e-10 && abs(ycp - Y(icp+1,jcp)) < 1e-10)
            l = indmapping(icp+1,jcp);
            Eh1(iter,l) = 1.0;
            Eh3(iter,l) = 1.0;
        elseif (abs(xcp - X(icp+1,jcp+1)) < 1e-10 && abs(ycp - Y(icp+1,jcp+1)) < 1e-10)
            l = indmapping(icp+1,jcp+1);
            Eh1(iter,l) = 1.0;
            Eh3(iter,l) = 1.0;
        else
        
            % Bi-linear interpolation matrix
            denom = -1/(xcp-X(icp,jcp)) + 1/(xcp-X(icp,jcp+1));
            wx = [-1/(xcp - X(icp,jcp)), 1/(xcp - X(icp,jcp+1))]/denom;
            denom = -1/(ycp-Y(icp,jcp)) + 1/(ycp-Y(icp+1,jcp));
            wy = [-1/(ycp - Y(icp,jcp)), 1/(ycp - Y(icp+1,jcp))]/denom;
            l = indmapping(icp,jcp);
            Eh1(iter,l) = wy(1)*wx(1);
            l = indmapping(icp,jcp+1);
            Eh1(iter,l) = wy(1)*wx(2);
            l = indmapping(icp+1,jcp);
            Eh1(iter,l) = wy(2)*wx(1);
            l = indmapping(icp+1,jcp+1);
            Eh1(iter,l) = wy(2)*wx(2);

            % Bi-cubic interpolation matrix
            denom = -1/3/(xcp-X(icp,jcp-1)) + 1/(xcp-X(icp,jcp)) + ...
                    -1/(xcp-X(icp,jcp+1)) + 1/3/(xcp-X(icp,jcp+2));
            wx = [-1/3/(xcp-X(icp,jcp-1)), 1/(xcp-X(icp,jcp)), ...
                  -1/(xcp-X(icp,jcp+1)), 1/3/(xcp-X(icp,jcp+2))]/denom;
            denom = -1/3/(ycp-Y(icp-1,jcp)) + 1/(ycp-Y(icp,jcp)) + ...
                    -1/(ycp-Y(icp+1,jcp)) + 1/3/(ycp-Y(icp+2,jcp));
            wy = [-1/3/(ycp-Y(icp-1,jcp)), 1/(ycp-Y(icp,jcp)), ...
                  -1/(ycp-Y(icp+1,jcp)), 1/3/(ycp-Y(icp+2,jcp))]/denom;
            for m=1:4
                for n=1:4
                    l = indmapping(icp-2+m,jcp-2+n);
                    Eh3(iter,l) = wy(m)*wx(n);
                end
            end
        end

        
        % Laplacian matrixes (note for edge nodes, row values dont
        % matter because Eh*Lh and Eh*v will not be affected... 
        % so set to identity
        if (cval(k) == 2)
            Lh(iter,iter) = 1.0;
            Acacp(iter,:) = -Eh3(iter,:);
            Acacp(iter,iter) = Acacp(iter,iter) + 1.0;
            bcacp(iter) = 0.0;
        else
            Acacp(iter,iter) = 1.0 + (1+kappa(k)*phi(k))*(4.0 + ...
                 0.5*(kappa(i,j+1)*phi(i,j+1) + kappa(i,j)*phi(i,j)) +...
                 0.5*(kappa(i,j-1)*phi(i,j-1) + kappa(i,j)*phi(i,j)) +...
                 0.5*(kappa(i+1,j)*phi(i+1,j) + kappa(i,j)*phi(i,j)) +...
                 0.5*(kappa(i-1,j)*phi(i-1,j) + kappa(i,j)*phi(i,j)) )/dx^2;

            l = indmapping(i,j+1);
            Acacp(iter,l) = -(1+kappa(k)*phi(k))*(1.0 + ...
                 0.5*(kappa(i,j+1)*phi(i,j+1) + kappa(i,j)*phi(i,j)) )/dx^2;

            l = indmapping(i,j-1);
            Acacp(iter,l) = -(1+kappa(k)*phi(k))*(1.0 + ...
                     0.5*(kappa(i,j-1)*phi(i,j-1) + kappa(i,j)*phi(i,j)) )/dx^2;

            l = indmapping(i+1,j);
            Acacp(iter,l) = -(1+kappa(k)*phi(k))*(1.0 + ...
                0.5*(kappa(i+1,j)*phi(i+1,j) + kappa(i,j)*phi(i,j)) )/dx^2;

            l = indmapping(i-1,j);
            Acacp(iter,l) = -(1+kappa(k)*phi(k))*(1.0 + ...
                0.5*(kappa(i-1,j)*phi(i-1,j) + kappa(i,j)*phi(i,j)) )/dx^2;
            
            bcacp(iter) = f(X(k)-phi(k)*nx(k),Y(k)-phi(k)*ny(k));
                        
            Lh(iter,iter) = -4.0/dx^2;
            l = indmapping(i,j+1);
            Lh(iter,l) = 1.0/dx^2;
            l = indmapping(i,j-1);
            Lh(iter,l) = 1.0/dx^2;
            l = indmapping(i+1,j);
            Lh(iter,l) = 1.0/dx^2;
            l = indmapping(i-1,j);
            Lh(iter,l) = 1.0/dx^2;
        end        
        
        % RHS evaluated at CP
        bcp(iter) = f(X(k)-phi(k)*nx(k),Y(k)-phi(k)*ny(k));
    end
    
    Acacp = sparse(Acacp);
    Lh = sparse(Lh);
    Eh1 = sparse(Eh1);
    Eh3 = sparse(Eh3);
    
    M = Eh1*Lh - 4/dx^2*(speye(L) - Eh3);
    A = speye(L) - M;
    uCP = A\bcp;
    
    unew = Acacp\bcacp;

    xcp = X(ind) - phi(ind).*nx(ind);
    ycp = Y(ind) - phi(ind).*ny(ind);
    uex = sin(atan2(ycp,xcp)) + sin(12*atan2(ycp,xcp));
    e2 = [sqrt(sum((uCP-uex).^2)/length(uex)) sqrt(sum((unew-uex).^2)/length(uex))];
    eI = [max(abs(uCP-uex)) max(abs(unew-uex))];
    
    fprintf('L-2 error (CP vs CACP): %e vs %e\n', e2(1), e2(2));
    fprintf('L-inf error (CP vs CACP): %e vs %e\n', eI(1), eI(2));
    fprintf('Non-zeros (CP vs CACP): %d vs %d\n\n', nnz(A), nnz(Acacp));
    
    Ucp = zeros(size(X));
    Ucp(ind) = uCP;
    Ucacp = zeros(size(X));
    Ucacp(ind) = unew;
    Uex = zeros(size(X));
    Uex(ind) = uex;
        
end

function val = RHS_circle(x,y)

    theta = atan2(y,x);
    val = 2*sin(theta) + 145*sin(12*theta);
    
end

function val = RHS_clover(x,y)

    global alpha

    theta = atan2(y,x);
    u = sin(theta) + sin(12*theta);
    up = cos(theta) + 12*cos(12*theta);
    upp = -sin(theta) - 144*sin(12*theta);
    g = 1 + alpha*cos(4*theta - pi);
    gp = -4.0*alpha*sin(4*theta - pi);
    gpp = -16.0*alpha*cos(4*theta - pi);
    
    val = -upp/(gp^2 + g^2) + (gp*gpp + g*gp)*up/(gp^2 + g^2)^2 + u;
    
end