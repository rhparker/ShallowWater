% perform one iteration of RK4 method
% assume the spatial operator is independent of t
%   k is time step
%   u is current discrete solution vector
%   f is function discrete spatial operator

function uout = RK4(u, k, f)
	k1 = k * f(u);
	k2 = k * f( u + 0.5 * k1 );
	k3 = k * f( u + 0.5 * k2 );
	k4 = k * f( u + k3 );
	uout =  u + (1.0/6) * (k1 + 2 * k2 + 2 * k3 + k4);
end