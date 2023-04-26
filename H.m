function H_result = H(u, tau)
    % Get the size of the input vector u
    n = length(u);
    
    % Compute the matrix product u*uH (Hermitian transpose of u)
    uuH = u * u';
    
    % Create the identity matrix of size n x n
    I = eye(n);
    
    % Compute the transformation I - (1/tau)*uuH
    H_result = I - (1/tau) * uuH;
end