
function B = bidiag_francis_step(B)
    [m, n] = size(B);
    if m <= 2
        return
    end
    
    % Compute the initial values of T11, T21, and Tmm
    T11 = B(1, 1)^2;
    T21 = B(1, 1) * B(1, 2);
    Tmm = B(2, 3)^2 + B(3, 3)^2;

    % Compute the Givens rotation that annihilates the (1,2) and (2,1) elements
    G = Givens_rotation( [T11- Tmm
                                 T21      ]);
    display(B);
    % Apply the rotation to the first two columns of B
    B(1:2, 1:2) = B(1:2, 1:2) * G;
    display(B);

    % Compute the updated values of T11, T21, and Tmm
    for i = 2:min(m-2, n-1)
        % Compute the updated values of T11, T21, and Tmm
        T11 = B(1, 1);
        T21 = B(2, 1);
        Tmm = B(3,3);
        
        T11_Tmm= T11-Tmm;
        % Compute the Givens rotation that annihilates the (i,i+1) and (i+1,i) elements
        G = Givens_rotation( [T11_Tmm
                                 T21      ]);

        % Apply the rotation to the (i,i+1) and (i+1,i) elements
        B(i-1:i, i-1:i) = G' * B(i-1:i, i-1:i);
        
        % Compute the updated values of T11, T21, and Tmm
        T11 = B(1, 1)^2;
        T21 = B(1, 1) * B(1, 2);
        Tmm = B(2, 3)^2 + B(3, 3)^2;

            % Compute the Givens rotation that introduces a new bulge in the matrix
             G = Givens_rotation( [T11- Tmm
                                 T21      ]);

            % Apply the rotation to the (i+1,i) and (i+1,i+1) elements
            B(i+1:i+2, i:i+1) = B(i+1:i+2, i:i+1) * G;
            B(i:i+1, i+1:i+2) = G' * B(i:i+1, i+1:i+2);
            
        end
        end
