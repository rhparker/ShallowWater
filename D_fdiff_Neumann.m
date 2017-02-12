function [D, D2, D3, D4, D5] = D_fdiff_Neumann(N, h)

% d_x
D = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N); 
D = (D - D');

D(1,2) = 0; 
D(N,N-1) = 0;

D = D./(2 * h);

% d_xx
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],ones(N,1),N,N);
D2 = (D2 + D2');

D2(1,2) = 2; 
D2(N,N-1) = 2;

D2 = D2./h^2;

% d_xxx
D3 = -2 * sparse(1:N-1,[2:N],ones(N-1,1),N,N) + sparse(1:N-2,[3:N],ones(N-2,1),N,N);
D3 = (D3 - D3');

% Neumann boundary conditions (third derivative 0 at L and R)
D3(1,2)   = 0; D3(1,3)   = 0;
D3(N,N-1) = 0; D3(N,N-2) = 0;
% Neumann boundary conditions (inherited from D)
D3(2,2)= -1; D3(N-1,N-1) = 1;
% scale by 2h^3
D3 = D3./(2 * h^3);

% d_xxxx
D4 = sparse(1:N,[1:N],3 * ones(N,1),N,N) - sparse(1:N-1,[2:N],4 * ones(N-1,1),N,N) + sparse(1:N-2,[3:N],ones(N-2,1),N,N);
D4 = (D4 + D4');

D4(1,2) = -8; D4(1,3) = 2;
D4(2,2) = 7;
D4(N,N-1) = -8; D4(N,N-2) = 2;
D4(N-1,N-1) = 7;

% scale by h^4
D4 = D4./h^4;

% d_xxxxx
D5 = 5 * sparse(1:N-1,[2:N],ones(N-1,1),N,N) -4 * sparse(1:N-2,[3:N],ones(N-2,1),N,N) + sparse(1:N-3,[4:N],ones(N-3,1),N,N);
D5 = (D5 - D5');

D5(1,2) = 0; D5(1,3) = 0; D5(1,4) = 0;
D5(2,2) = 4;  D5(2,3) = 4;
D5(3,2) = -6;
D5(N,N-1) = 0; D5(N,N-2) = 0; D5(N,N-3) = 0;
D5(N-1, N-1) = -4; D5(N-1, N-2) = -4;
D5(N-2, N-1) = 6;

% scale by h^5
D5 = D5./(2 * h^5);

end