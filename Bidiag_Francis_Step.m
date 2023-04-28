% function Bi_next = Bidiag_Francis_Step(B)
% 
% m=size(B,1);
% T11=B(1,1);
% display(T11);
% T11=(T11*T11);
% 
% T21=B(1,2)*B(1,1);
% display(T21);
% Tmm=(B(m-1,m)^2)+ (B(m,m)^2);
% display(Tmm);
% 
% for i=1,m
% tau = B(i,i)^2 + B(i+1,i)^2;
% sigma = B(i+1,i)/tau;
% gamma = B(i,i)/tau;
% 
% G0=[gamma ,sigma;-sigma gamma];
% 
% Bvector=[B(i,i);B(i+1,i)];
% 
% Bvector=G0'*Bvector;
% display(B);
% display(Bvector);
% 
% 
% end


function Bi_next = Bidiag_Francis_Step(B)
m = size(B,1);
T11 = B(1,1)^2;
T21 = B(1,2)*B(1,1);
Tmm = B(m-1,m)*B(m,m) + B(m,m)^2;

% Compute first Givens rotation
tau = T11 + T21;
sigma = B(1,2)/tau;
gamma = B(1,1)/tau;
G = [gamma, sigma; -sigma, gamma];

% Apply Givens rotation to introduce the bulge
B(1:2, 1:2) = G * B(1:2, 1:2);

% Apply Givens rotations to chase the bulge to the right
for i = 1:m-2
    % Compute Givens rotation parameters
    tau = B(i,i)^2 + B(i+1,i)^2;
    sigma = B(i+1,i)/tau;
    gamma = B(i,i)/tau;
    
    % Compute and apply left Givens rotation
    G = [gamma, sigma; -sigma, gamma];
    B(i:i+1, i:m) = G' * B(i:i+1, i:m);
    
    % Compute and apply right Givens rotation
    tau = B(i,i+1)^2 + B(i,i+2)^2;
    sigma = B(i,i+2)/tau;
    gamma = B(i,i+1)/tau;
    G = [gamma, sigma; -sigma, gamma];
    B(1:i+2, i+1:i+2) = B(1:i+2, i+1:i+2) * G;
end

% Apply final Givens rotation to chase the bulge off the matrix
tau = B(m-1,m-1)^2 + Tmm;
sigma = B(m-1,m)/tau;
gamma = B(m-1,m-1)/tau;
G = [gamma, sigma; -sigma, gamma];
B(m-1:m, :) = G' * B(m-1:m, :);

Bi_next = B;
end
