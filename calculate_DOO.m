% SOM is the stripped observability matrix
function D = calculate_DOO(SOM)
% Set λ as a small positive real number, which deals with the situation that R is rank deficient and R'R is therefore with unstable matrix inversion
lamda = 1e-10; 
% QR decomposition
R = triu(qr(SOM)); 
% Calculate O+O
OO = (R'*R + lamda*eye(size(SOM,2)))\(R'*R);
% Get the DOO of each state
D = sum(OO.*OO,2);
end
