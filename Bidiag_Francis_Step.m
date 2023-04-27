function Bi = bidiag_francis_step(Bi)
% Perform a single step of the Bidiagonal Francis Step algorithm on the
% bidiagonal matrix Bi, working on the i-th superdiagonal element.
% Input:
%   - Bi: a bidiagonal matrix of size m x n (m >= n)
% Output:
%   - Bi: the updated bidiagonal matrix

% Determine the size of the matrix Bi
[m, n] = size(Bi);

% Loop over the superdiagonal elements of Bi
for i = 1:min(m-1, n)
    % Extract the i-th superdiagonal element and the (i+1)-th element in the
    % same column
    a = Bi(i, i);
    c = Bi(i+1, i);

    % Compute the rotation matrix that annihilates the two elements
    x = sqrt(a^2 + c^2);
    cos_theta = a / x;
    sin_theta = -c / x;
    G = [cos_theta, -sin_theta; sin_theta, cos_theta];

    % Apply the rotation to the i-th and (i+1)-th rows of Bi
    Bi(i:i+1, :) = G * Bi(i:i+1, :);

    % Apply the transpose of the rotation to the i-th and (i+1)-th columns of Bi
    Bi(:, i:i+1) = Bi(:, i:i+1) * G';

    % Update the (i+1)-th superdiagonal element, if it exists
    if i < n-1
        % Extract the (i+1)-th superdiagonal element and the (i+2)-th element in the
        % same column
        b = Bi(i+1, i+1);
        d = Bi(i+1, i+2);

        % Compute the rotation matrix that introduces a new bulge in the matrix
        y = sqrt(b^2 + d^2);
        cos_phi = b / y;
        sin_phi = -d / y;
        H = [cos_phi, -sin_phi; sin_phi, cos_phi];

        % Apply the rotation to the (i+1)-th and (i+2)-th columns of Bi
        Bi(:, i+1:i+2) = Bi(:, i+1:i+2) * H;

        % Apply the transpose of the rotation to the (i+1)-th and (i+2)-th rows of Bi
        Bi(i+1:i+2, :) = H' * Bi(i+1:i+2, :);
    end
end
