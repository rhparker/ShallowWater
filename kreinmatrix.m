% computes the Krein matrix

% standard double pulses
load 100F_double2;
N = length(xout);

% % large domain for crossing
% load 1500F_domain_105;
% index = 4;
% xout = output_x(:,index);
% uout = output_u(:,index);
% par.c = output_u(end,index);
% L = output_L(index);
% N = length(xout);

% % differentiation operators
[D, ~, ~, ~, ~] = D_fourier(N, L);
% Jacobian
config.symmetry = 'none';
J_int = get_jacobian(xout, uout(1:end-1), par, config, 'integrated', 0);

% load 1500F_double2_1195;

% we need a basis for our subspace which extended 
% the kernel of the two operators
% (these are already orthogonal)

u = uout(1:end-1);
Du = D*u;

% pulse u and its derivative
v1 = u/norm(u);
v2 = Du/norm(Du);

% constant function, should be kernel of diff operator
v3 = ones(N,1) / norm(ones(N,1));

% turns out that ker D has two things in it, likely due to
% something about fourier spectral differentiation; in any
% case will use that as kernel
kerD = null(D);

% easiest way to get a basis is to tack on standard
% basis vectors and use Gram-Schmidt. QR accomplishes
% this, interestingly enough; the columns of M are
% an ON basis extending v1, etc 
Id = eye(N);

% % this basis extends v2 = Du
% M1 = [v2 Id(:,1:end-1)];
% [M1, R] = qr(M1, 0);
% 
% % this basis extends {v1, ker D}
% M2 = [v1 kerD Id(:,3:end-1)];
% [M2, R] = qr(M2, 0);

% this basis extends {v1, v2, Ker D}
M = [v1 v2 kerD Id(:,4:end-1)];
[M, ~] = qr(M, 0);

% projection matrix in M-coordinates
% projection gets rid of v1, v2, Ker D
% PM is proj matrix in M-coord, P in std coord
PM = Id; 
PM(:,1:4) = zeros(N, 4);
P = M*PM*M';

Hplus = J_int;
Hminus = -D*J_int*D;

% now we define R and Sinverse, using the formulas
%    R    = PM Hplus PM
%    Sinv = PM Hminus PM

R = P*Hplus*P;
Sinv = P*Hminus*P;

% we need to invert S, but only on the N-4 dimensional
% subspace which doesn't contain its kernel

Z = M'*Sinv*M;
Z2 = inv( Z(5:end,5:end) );
Z(5:end,5:end) = Z2;
S = P*M*Z*M'*P;

% another version of S, this time it acts as the identity on 
% the 4-dimensional kernel of the projection P; this is
% more convenient for finding negative eigenvalues
% since we don't have any eigenvalues of zero; we checked
% to make sure the negative eigenvalues are the same
% as those of S.
Z = M'*Sinv*M;
Z(:,1:4) = Id(:,1:4);
Snonsingular = inv( M*Z*M' );

% eigenvalues and eigenfunctions of S (using Snonsingular)
[V, lambda] = eig(Snonsingular);
lambda = diag(lambda);
% find negative eigenvalues and corresp eigenfunctions
lNeg = lambda(find(lambda < 0));
vNeg = V(:,find(lambda < 0));
nNeg = length(lNeg);

% plot the negative eigenfunctions
figure;
plot(xout, vNeg);
title('Eigenfunctions corresponding to negative eigenvalues of the operator S');
legendCell = cellstr(num2str(lNeg, 'nu=%-d'))
legend(legendCell);

% plot derivatives of these
dvNeg = D*vNeg;
for index = 1:length(lNeg)
    dvNeg(:,index) = dvNeg(:,index) / norm(dvNeg(:,index));
end
figure;
plot(xout, dvNeg);
title('Derivatives of eigenfunctions corresponding to negative eigenvalues of the operator S');
legendCell = cellstr(num2str(lNeg, 'nu=%-d'))
legend(legendCell);

% now we will need projection on orthogonal complement
% of negative eigenspace of S

% this basis extends {vNeg}
MNeg = [vNeg Id(:,nNeg:end-1)];
[MNeg, ~] = qr(MNeg, 0);

% this basis extends {v1, v2, kerD, vNeg}
% will need it to construct the resolvent operator
% i.e. invert (R2 - z S2)
M2 = [v1 v2 kerD vNeg Id(:,nNeg + 5:end)];
[M2, ~] = qr(M2, 0);

% projection gets rid of vNeg
PNegPerp = Id; 
PNegPerp(:,1:nNeg) = zeros(N, nNeg);
PNegPerp = MNeg*PNegPerp*MNeg';

% define R2 and S2 as projections on perp space to vNeg
R2 = PNegPerp*R*PNegPerp;
S2 = PNegPerp*S*PNegPerp;

% make the Rhat matrix, expect to be diagonal since
% the two functions in vNeg are odd and even
Rhat = zeros(nNeg, nNeg);
for i = 1:nNeg
    for j = 1:nNeg
        Rhat(i,j) = vNeg(:,i)' * (R*vNeg(:,j));
    end
end

% negative eigenvalues on diagonal of matrix
ND = diag(lNeg);

% from this point on, things will depend on z

z = 0.0691^2;

% % see if we care about Rtilde (answer = no)
% Sroot = sqrtm(S2);
% Rtilde = inv(Sroot)*R2*inv(Sroot);
% resolv2 = inv(Sroot)*inv(Rtilde - z*Id)*inv(Sroot);

% compute resolvent operator
Z = M2'*(R2 - z*S2)*M2;
Z2 = inv( Z(5 + nNeg:end, 5+ nNeg:end) );
Z(5 + nNeg:end,5 + nNeg:end) = Z2;
resolv = M2*Z*M2';

% Cz portion of Krein matrix
Cz = zeros(nNeg, nNeg);
for i = 1:nNeg
    for j = 1:nNeg
        Cz(i,j) = (PNegPerp*resolv*PNegPerp*R*vNeg(:,i))' * (PNegPerp*R*vNeg(:,j));
    end
end

Kz = Rhat - z*ND - Cz;

% save things so we can compute this easily for arbitrary z
% save 1500F_double2_1195_krein M2 R2 S2 R S lNeg vNeg PNegPerp Rhat ND






