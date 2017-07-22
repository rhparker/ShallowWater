% computes the Krein matrix, or at least tries to

load 100F_double2;

N = length(xout);

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

% turns out that D has two things in it, likely due to
% something about fourier spectral differentiation
kerD = null(D);

% easiest way to get a basis is to tack on standard
% basis vectors and use Gram-Schmidt. QR accomplishes
% this, interestingly enough; the columns of M are
% an ON basis extending v1, etc 
Id = eye(N);

% this basis extends v2 = Du
M1 = [v2 Id(:,1:end-1)];
[M1, R] = qr(M1, 0);

% this basis extends {v1, ker D}
M2 = [v1 kerD Id(:,3:end-1)];
[M2, R] = qr(M2, 0);

% this basis extends {v1, v2, Ker D}
M = [v1 v2 kerD Id(:,4:end-1)];
[M, R] = qr(M, 0);

% % projection matrix in M-coordinates
% % projection gets rid of v1, v2, Ker D
PM = Id; 
PM(:,1:4) = zeros(N, 4);
P = M*PM*M';

Hplus = J_int;
Hminus = -D*J_int*D;

% We want these to be invertible, so in theory we 
% restrict to perpendicular space to kernel of Hplus/Hminus
% then invert. Practically, this is annoying, so 
% we can get rid of the kernel my modifying
% Hminus to act as the identity on its kernel

% for now, we get rid of a 4-dim subspace

% do this to Hplus to get R
Z = M'*Hplus*M;
Z(:,1:4) = Id(:,1:4);
R = M*Z*M';

% do this to Hminus to get SInv, invert to get S
Z = M'*Hminus*M;
Z(:,1:4) = Id(:,1:4);
SInv = M*Z*M';
S = inv(SInv);

% eigenvalues and eigenfunctions of S
[V, lambda] = eig(S);
lambda = diag(lambda);
% find negative eigenvalues and corresp eigenfunctions
lNeg = lambda(find(lambda < 0));
vNeg = V(:,find(lambda < 0));

% now we will need projection on orthogonal complement
% of negative eigenspace of S

% this basis extends {vNeg}
MNeg = [vNeg Id(:,2:end-1)];
[MNeg, R] = qr(MNeg, 0);

% projection gets rid of vNeg
PNegPerp = Id; 
PNegPerp(:,1:2) = zeros(N, 2);
PNegPerp = M*PNegPerp*M';

% we also want a version of R2 = PNeg R PNeg
% as before, since we don't want to reduce the rank,
% we will make R2 act as the identity on the vNeg

Z = MNeg'*R*MNeg;
Z(:,1:2) = Id(:,1:2);
R2 = MNeg*Z*MNeg';

% do the same thing for S2 = PNeg S PNeg

Z = MNeg'*S*MNeg;
Z(:,1:2) = Id(:,1:2);
S2 = MNeg*Z*MNeg';

% make the Rhat matrix, expect to be diagonal since
% the two functions in vNeg are odd and even
Rhat = zeros(2, 2);
for i = 1:2
    for j = 1:2
%         Rhat(i,j) = vNeg(:,i)' * (R*vNeg(:,j));
        Rhat(i,j) = vNeg(:,i)'* (PM*R*PM*vNeg(:,j));
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

resolv = inv(R2 - z*S2);

Cz = zeros(2, 2);
for i = 1:2
    for j = 1:2
        Cz(i,j) = (PNegPerp*resolv*PNegPerp*R*vNeg(:,i))' * (PNegPerp*R*vNeg(:,j));
    end
end

Kz = Rhat - z*ND - Cz;






