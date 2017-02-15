function [F,J] = FixedPointSecantPredictorCorrector(v,v1,v0,N,config,D,D2,D3,D4,D5,par,contPar)

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

F = zeros(size(v));
[Fu,Ju] = integratedequation(u,par,N,config,D,D2,D3,D4,D5);          % evaluate system at predictor u and p

F(1:end-1) = Fu;
F(end)    = alpha' * (u - u1) + beta * (p - p1) - contPar.ds; % 

if nargout > 1
    h = 1e-8;
    parpe = setfield(par,contPar.Name,p+h);
    
    Fu2 = integratedequation(u,parpe,N,config,D,D2,D3,D4,D5);
    
    Fp = (Fu2 - Fu)/h;
    
    J = [Ju     Fp;
         alpha' beta];
     
    J = sparse(J);
end

