function [e2,eI,numnz,condn,X,Y,Z,Ucp,Ucacp,Uex] = shiftp_equation3(M, Nthreads) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function to solve an elliptic PDE on a sphere in R^3 using the Closest 
%   Point and Curvature-Augmented Closest Point methods.
%
%   M           - number of grid cells in one direction
%   Nthreads    - number of parallel workers to use in building the
%                  matrices
%
%   e2          - L2 error for both methods
%   eI          - L-inf error for both methods
%   numnz       - number of non-zeros for both methods
%   condn       - condition number estimate for both methods
%   (X,Y,Z)     - meshgrid for plotting purposes
%   Ucp         - Closest Point method solution
%   Ucacp       - Curvature-Augmented Closest Point solution
%   Uex         - Exact solution
%
%   See README.md for link to paper with more information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    f = @(x,y,z)(RHS_sphere(x,y,z));
    sdf = @(x,y,z)(sdf_sphere(x,y,z));
    
    %initialize grid
    [X,Y,Z] = meshgrid(linspace(-2,2,M+1),linspace(-2,2,M+1),linspace(-2,2,M+1));
    dx = X(1,2,1)-X(1,1,1);
    fprintf('\ndx = %f\n',dx);

    [Phi, Nx, Ny, Nz, Kappa, H] = sdf(X,Y,Z);
     
    % Classify interface grid points
    cval = classify(Phi);
    ind = find(cval > 0);
    L = length(ind);
    
    % Generate mapping from (i,j) linear system index & convert sdf
    % variables
    indmapping = zeros(size(X));
    phi = zeros(L,1);
    nx = zeros(L,1);
    ny = zeros(L,1);
    nz = zeros(L,1);
    kappa = zeros(L,1);
    h = zeros(L,1);
    for iter = 1:L
        [i,j,k] = ind2sub(size(X),ind(iter));
        indmapping(i,j,k) = iter;
        phi(iter) = Phi(i,j,k);
        nx(iter) = Nx(i,j,k);
        ny(iter) = Ny(i,j,k);
        nz(iter) = Nz(i,j,k);
        kappa(iter) = Kappa(i,j,k);
        h(iter) = H(i,j,k);
    end
    
    % Generate closest point information and delete unneeded variables for
    % memory
    xcp = X(ind) - Phi(ind).*Nx(ind);
    ycp = Y(ind) - Phi(ind).*Ny(ind);
    zcp = Z(ind) - Phi(ind).*Nz(ind);
    clear Phi Nx Ny Nz Kappa H
    
    % Parallelize if specified
    if (nargin >= 2)
        Lbreak = [0 (1:Nthreads-1)*floor(L/Nthreads) L];
    
        pool = gcp('nocreate');
        if (isempty(pool))
            pool = parpool(Nthreads);
        else
            delete(pool);
            pool = parpool(Nthreads);
        end

        fprintf('Building matrices using %d threads with %d points each...\n',Nthreads,floor(L/Nthreads));
        tic;
        for iter=1:Nthreads
            F(iter) = parfeval(pool,@buildSystem,6,Lbreak(iter)+1,Lbreak(iter+1),X(1,:,1),Y(:,1,1),Z(1,1,:),cval,ind,indmapping,phi,nx,ny,nz,kappa,h,f);
        end

        [Eh1,Eh3,Lh,Acacp,bcp,bcacp] = fetchOutputs(F(1));

        for iter=2:Nthreads
            [Eh1add,Eh3add,Lhadd,Ashadd,badd,bsadd] = fetchOutputs(F(iter));
            Eh1 = [Eh1; Eh1add];
            Eh3 = [Eh3; Eh3add];
            Lh = [Lh; Lhadd];
            Acacp = [Acacp; Ashadd];
            bcp = [bcp; badd];
            bcacp = [bcacp; bsadd];
        end
        toc;
    else
        fprintf('Building matrices using 1 thread...\n');
        tic;
        [Eh1,Eh3,Lh,Acacp,bcp,bcacp] = buildSystem(1,L,X(1,:,1),Y(:,1,1),Z(1,1,:),cval,ind,indmapping,phi,nx,ny,nz,kappa,h,f);
        toc;
    end

    Acacp = sparse(Acacp(:,1),Acacp(:,2),Acacp(:,3),L,L);
    Lh = sparse(Lh(:,1),Lh(:,2),Lh(:,3),L,L);
    Eh1 = sparse(Eh1(:,1),Eh1(:,2),Eh1(:,3),L,L);
    Eh3 = sparse(Eh3(:,1),Eh3(:,2),Eh3(:,3),L,L);

    fprintf('Finished building matrices, now solving...\n');
    M = Eh1*Lh - 6/dx^2*(speye(L) - Eh3);
    A = speye(L) - M;
    ucp = A\bcp;

    ucacp = Acacp\bcacp;
            
    theta = atan2(ycp,xcp);
    psi = atan2(sqrt(xcp.^2 + ycp.^2),zcp);
    uex = cos(3*theta).*sin(psi).^3.*(9*cos(psi).^2 - 1);    
    e2 = [sqrt(sum((ucp-uex).^2)/length(uex)) sqrt(sum((ucacp-uex).^2)/length(uex))];
    eI = [max(abs(ucp-uex)) max(abs(ucacp-uex))];
    numnz = [nnz(A) nnz(Acacp)];
    condn = [condest(A) condest(Acacp)];

    Ucp = zeros(size(X));
    Ucp(ind) = ucp;
    Ucacp = zeros(size(X));
    Ucacp(ind) = ucacp;
    Uex = zeros(size(X));
    Uex(ind) = uex;
        
end

function val = RHS_sphere(x,y,z)

    theta = atan2(y,x);
    psi = atan2(sqrt(x.^2 + y.^2),z);
    val = 31*cos(3*theta).*sin(psi).^3.*(9*cos(psi).^2 - 1);
    
end

function [Eh1,Eh3,Lh,Acacp,bcp,bcacp] = buildSystem(start,stop,x,y,z,cval,ind,indmapping,phi,nx,ny,nz,kappa,h,f)

    dx = x(2)-z(1);
    L = stop-start-1;
    
    Eh1 = zeros(L*8,3);
    nEh1 = 0;
    Eh3 = zeros(L*64,3);
    nEh3 = 0;
    Lh = zeros(L*7,3);
    nLh = 0;
    Acacp = zeros(L*7,3);
    nAsh = 0;
    bcp = zeros(L,1);
    bcacp = zeros(L,1);
    for iter=start:stop
        [i,j,k] = ind2sub([length(x),length(y),length(z)],ind(iter));

        % Interpolation matrices (need eps to avoid NaN in interpolants)
        xcp = x(j) - phi(iter)*nx(iter) + eps;
        ycp = y(i) - phi(iter)*ny(iter) + eps;
        zcp = z(k) - phi(iter)*nz(iter) + eps;
        icp = find(y > ycp,1) - 1;
        jcp = find(x > xcp,1) - 1;
        kcp = find(z > zcp,1) - 1;
        
        on_interface = 0;
        for s1 = [0 1]
            for s2 = [0 1]
                for s3 = [0 1]
                    indadj = [s1, s2, s3];
                    if (abs(xcp - x(jcp+indadj(2))) < 1e-10 && ...
                        abs(ycp - y(icp+indadj(1))) < 1e-10  && ...
                        abs(zcp - z(kcp+indadj(3))) < 1e-10)
                        l = indmapping(icp+indadj(1),jcp+indadj(2),kcp+indadj(3));
                        nEh1 = nEh1 + 1;
                        Eh1(nEh1,:) = [iter,l,1.0];
                        nEh3 = nEh3 + 1;
                        Eh3(nEh3,:) = [iter,l,1.0];
                        on_interface = 1;
                    end
                end
            end
        end
        
        if (on_interface == 0)
            % Bi-linear interpolation matrix
            denom = -1/(xcp-x(jcp)) + 1/(xcp-x(jcp+1));
            wx = [-1/(xcp - x(jcp)), 1/(xcp - x(jcp+1))]/denom;
            denom = -1/(ycp-y(icp)) + 1/(ycp-y(icp+1));
            wy = [-1/(ycp - y(icp)), 1/(ycp - y(icp+1))]/denom;
            denom = -1/(zcp-z(kcp)) + 1/(zcp-z(kcp+1));
            wz = [-1/(zcp - z(kcp)), 1/(zcp - z(kcp+1))]/denom;

            for m=1:2
                for n=1:2
                    for o=1:2
                        indadj = [m-1, n-1, o-1];
                        l = indmapping(icp+indadj(1),jcp+indadj(2),kcp+indadj(3));
                        nEh1 = nEh1 + 1;
                        Eh1(nEh1,:) = [iter,l,wy(m)*wx(n)*wz(o)];
                    end
                end
            end

            % Bi-cubic interpolation matrix
            denom = -1/3/(xcp-x(jcp-1)) + 1/(xcp-x(jcp)) + ...
                    -1/(xcp-x(jcp+1)) + 1/3/(xcp-x(jcp+2));
            wx = [-1/3/(xcp-x(jcp-1)), 1/(xcp-x(jcp)), ...
                  -1/(xcp-x(jcp+1)), 1/3/(xcp-x(jcp+2))]/denom;
            denom = -1/3/(ycp-y(icp-1)) + 1/(ycp-y(icp)) + ...
                    -1/(ycp-y(icp+1)) + 1/3/(ycp-y(icp+2));
            wy = [-1/3/(ycp-y(icp-1)), 1/(ycp-y(icp)), ...
                  -1/(ycp-y(icp+1)), 1/3/(ycp-y(icp+2))]/denom;
            denom = -1/3/(zcp-z(kcp-1)) + 1/(zcp-z(kcp)) + ...
                    -1/(zcp-z(kcp+1)) + 1/3/(zcp-z(kcp+2));
            wz = [-1/3/(zcp-z(kcp-1)), 1/(zcp-z(kcp)), ...
                  -1/(zcp-z(kcp+1)), 1/3/(zcp-z(kcp+2))]/denom;
            for m=1:4
                for n=1:4
                    for o=1:4
                        indadj = [m-2, n-2, o-2];
                        l = indmapping(icp+indadj(1),jcp+indadj(2),kcp+indadj(3));
                        nEh3 = nEh3 + 1;
                        Eh3(nEh3,:) = [iter,l,wy(m)*wx(n)*wz(o)];
                    end
                end
            end
        end

        
        % Laplacian matrixes (note for edge nodes, row values dont
        % matter because Eh*Lh and Eh*v will not be affected... 
        % so set to identity
        if (cval(i,j,k) == 2)
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,iter,1.0];
            tmpind = nEh3;
            flag = 0;
            while (tmpind > 0 && Eh3(tmpind,1) == iter)
                nAsh = nAsh + 1;
                Acacp(nAsh,:) = [Eh3(tmpind,1:2), -6/dx^2*Eh3(tmpind,3)];
                if (Eh3(tmpind,2) == iter)
                    Acacp(nAsh,3) = Acacp(nAsh,3) + 6/dx^2;
                    flag = 1;
                end
                tmpind = tmpind - 1;
            end
            if (flag == 0)
                nAsh = nAsh + 1;
                Acacp(nAsh,:) = [iter, iter, 6/dx^2];
            end
            bcacp(iter-start+1) = 0.0;
        else
            m = indmapping(i,j,k);
            n = indmapping(i-1,j,k);
            simh = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            n = indmapping(i+1,j,k);        
            siph = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            n = indmapping(i,j-1,k);
            sjmh = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            n = indmapping(i,j+1,k);
            sjph = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            n = indmapping(i,j,k-1);
            skmh = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            n = indmapping(i,j,k+1);
            skph = 0.5*((1.0 + kappa(n)*phi(n))/(1.0 + h(n)*phi(n)) + ...
                        (1.0 + kappa(m)*phi(m))/(1.0 + h(m)*phi(m)));
            s = (1.0 + kappa(m)*phi(m))*(1.0 + h(m)*phi(m));

            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,iter,1.0+s*(simh + siph + sjmh + sjph + skmh + skph)/dx^2];

            l = indmapping(i,j+1,k);
            nAsh = nAsh+1;
            Acacp(nAsh,:) = [iter,l,-s*sjph/dx^2];

            l = indmapping(i,j-1,k);
            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,l,-s*sjmh/dx^2];

            l = indmapping(i+1,j,k);
            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,l,-s*siph/dx^2];

            l = indmapping(i-1,j,k);
            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,l,-s*simh/dx^2];

            l = indmapping(i,j,k+1);
            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,l,-s*skph/dx^2];

            l = indmapping(i,j,k-1);
            nAsh = nAsh + 1;
            Acacp(nAsh,:) = [iter,l,-s*skmh/dx^2];
            
            bcacp(iter-start+1) = f(xcp,ycp,zcp);    
            
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,iter,-6.0/dx^2];
            l = indmapping(i,j+1,k);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
            l = indmapping(i,j-1,k);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
            l = indmapping(i+1,j,k);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
            l = indmapping(i-1,j,k);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
            l = indmapping(i,j,k+1);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
            l = indmapping(i,j,k-1);
            nLh = nLh + 1;
            Lh(nLh,:) = [iter,l,1.0/dx^2];
        end
        
        % RHS evaluated at CP
        bcp(iter-start+1) = f(xcp,ycp,zcp);
    end        
    
    Eh1 = Eh1(1:nEh1,:);
    Eh3 = Eh3(1:nEh3,:);
    Lh = Lh(1:nLh,:);
    Acacp = Acacp(1:nAsh,:);
end
