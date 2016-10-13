function [val] = classify(phi)

    if(length(size(phi)) == 2)
        [val] = classify2D(phi);
    elseif (length(size(phi)) == 3)
        [val] = classify3D(phi);
    else
        error('Problem classifying grid nodes');
    end

end

function [val] = classify2D(phi)

    % Mark with neighbors
    interface = phi(2:end-1,2:end-1).*phi(2:end-1,3:end) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(2:end-1,1:end-2) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(3:end,2:end-1) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(1:end-2,2:end-1) <= eps;
    
    % Mark with corners
    interface = interface | phi(2:end-1,2:end-1).*phi(3:end,3:end) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(3:end,1:end-2) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(1:end-2,3:end) <= eps;
    interface = interface | phi(2:end-1,2:end-1).*phi(1:end-2,1:end-2) <= eps;
    
    %%%% Classify grid nodes %%%%
    % value of 1 indicates interpolation node
    % value of 2 indicates edge node
    val = interface;
    val = padarray(val,[1 1]);
    val = val + 0; % convert to double from logical
    
    ind = find(val == 1);
    for iter = 1:length(ind)
        k = ind(iter);
        [i,j] = ind2sub(size(phi),k);
        if (val(i,j+1) == 0)
            val(i,j+1) = 1;
        end
        if (val(i,j-1) == 0)
            val(i,j-1) = 1;
        end
        if (val(i+1,j) == 0)
            val(i+1,j) = 1;
        end
        if (val(i-1,j) == 0)
            val(i-1,j) = 1;
        end
        if (val(i+1,j+1) == 0)
            val(i+1,j+1) = 1;
        end
        if (val(i-1,j+1) == 0)
            val(i-1,j+1) = 1;
        end
        if (val(i+1,j-1) == 0)
            val(i+1,j-1) = 1;
        end
        if (val(i-1,j-1) == 0)
            val(i-1,j-1) = 1;
        end
    end
    
    % Classify edge nodes
    ind = find(val == 1);
    for iter=1:length(ind)
        k = ind(iter);
        [i,j] = ind2sub(size(phi),k);
        
        if (val(i,j+1) == 0) 
            val(i,j+1) = 2;
        end
        if (val(i,j-1) == 0) 
            val(i,j-1) = 2;
        end
        if (val(i+1,j) == 0) 
            val(i+1,j) = 2;
        end
        if (val(i-1,j) == 0) 
            val(i-1,j) = 2;
        end
    end
end

function [val] = classify3D(phi)

    interface = isnan(phi(2:end-1,2:end-1,2:end-1));

    % Mark interface cells
    for s1 = -1:1
        for s2 = -1:1
            for s3 = -1:1
                indadj = [s1, s2, s3];
                interface = interface | ...
                phi(2:end-1,2:end-1,2:end-1).*phi((2:end-1) + indadj(1), (2:end-1) + indadj(2), (2:end-1) + indadj(3)) <= 0.0;        
            end
       end
    end
        
    %%%% Classify grid nodes %%%%
    % value of 1 indicates interpolation node
    % value of 2 indicates edge node
    val = interface;
    val = padarray(val,[1 1 1]);
    val = val + 0; % convert to double from logical
    
    % Mark interpolation nodes
    ind = find(val == 1);
    for iter = 1:length(ind)
        l = ind(iter);
        [i,j,k] = ind2sub(size(phi),l);
        
        for s1 = -1:1
            for s2 = -1:1
                for s3 = -1:1
                    indadj = [s1, s2, s3];
                    if (val(i+indadj(1),j+indadj(2),k+indadj(3)) == 0)
                        val(i+indadj(1),j+indadj(2),k+indadj(3)) = 1;
                    end
                end
           end
        end        
    end
    
    % Classify edge nodes
    ind = find(val == 1);
    for iter=1:length(ind)
        l = ind(iter);
        [i,j,k] = ind2sub(size(phi),l);
        
        for off = 1:3
            for s = [-1 1]
                indadj = zeros(1,3);
                indadj(off) = s;
                if (val(i+indadj(1),j+indadj(2),k+indadj(3)) == 0)
                    val(i+indadj(1),j+indadj(2),k+indadj(3)) = 2;
                end
            end
        end
    end
end


