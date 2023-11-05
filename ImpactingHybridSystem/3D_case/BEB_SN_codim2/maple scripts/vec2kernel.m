%> the aim here is to calculate the tangent space for given vector
%> the numerical methods to be used here is Gram-Schmit Orthogonalization

function Kel = vec2kernel(v) 
    v   = v/norm(v);
    dim = length(v);
    
    
    %> kernel
    Kel = [];
    Kel = [Kel;v];
    %> use the Gram-Schimit Orthogobalization process
    for  i = 2 : dim
        %> initialize a vector with random vector which may composed of the
        %> full basis
        tmp     = rand(dim,1);
        tmp     = tmp./norm(tmp);
        for j = 1: i-1
            v_tmp   = Kel(:,j);
            tmp     = proj_deduct(tmp,v_tmp);
        end
         Kel = [Kel, tmp/norm(tmp)];

    end
    
    Kel = Kel(:,2:end);

    %> FUNCTION
    function tmp = proj_deduct(tmp,v)
        %> tmp is the temporary vector 
        %> v is the vector to be deducted 
            tmp = tmp - (tmp'*v)/norm(v)/norm(v)*v;
    end
   
end