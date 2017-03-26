function [F,J] = FixedPointSecantPredictorCorrector(x,v,v1,v0,eq,N,config,D,D2,D3,D4,D5,par,contPar,vsymm)

u = v(1:end-1); % Predictor
p = v(end); 

u1 = v1(1:end-1); % second point
p1 = v1(end); 

u0 = v0(1:end-1); % first point
p0 = v0(end);

d = norm( v1 - v0 , 2 );

alpha = (u1 - u0) / d;
beta  = (p1 - p0) / d;

par = setfield(par,contPar.Name,p); % update par to predictor p

[Fu,Ju] = eq(x,u,par,N,config,D,D2,D3,D4,D5,vsymm);          % evaluate system at predictor u and p

% F = zeros(size(v));
% F(1:end-1) = Fu;
% F(end)    = alpha' * (u - u1) + beta * (p - p1) - contPar.ds; 

F = [Fu; alpha' * (u - u1) + beta * (p - p1) - contPar.ds]; 


if nargout > 1
    h = 1e-8;
    parpe = setfield(par,contPar.Name,p+h);
    
    Fu2 = eq(x,u,parpe,N,config,D,D2,D3,D4,D5,vsymm);
    
    Fp = (Fu2 - Fu)/h;
    
    J = [Ju     Fp;
         alpha' beta];
     
    J = sparse(J);
end

