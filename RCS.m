function y = RCS(dist_fac,theta,phi,N,lambda,epsi,a) %[Ex,Ey,Ez] 

% Note: far field distance r>> sphere radius a. Ideally, choose r >= 10*a, so dist_fac>=10
% theta[elevation] can be varied from (0,pi) and phi[azimuth] from (0,2pi)
% N: upper cutoff for summation in Mie series scattered field
% lambda: Wavelength in metres of incident wave
% epsi: Complex permittivity of sphere

k = (2*pi)/lambda;  	  % free space wavevector 

r = dist_fac*a;		  % Far field radial distance r = distance factor*sphere radius

m = sqrt(epsi);      % refractive index

rho = k*r;          % scaled far field distance
x = k*a;            % scaled sphere radius

sumr = 0; sumt = 0; sump = 0; % r,theta,phi components of E_scat initialized to 0

% Pn1(cos t)/sin t ---> pi_n 
% dPn1(cos t)/dt -----> tau_n

mu = cos(theta);
pai = zeros(1,N); tau = zeros(1,N);

pai(1) = 1; tau(1) = mu;
pai(2) = 3*mu; tau(2) = 6*mu^2 - 3;

for n = 3:N
    pai(n) = ( (2*n-1)*mu*pai(n-1) - n*pai(n-2) )/(n-1);
    tau(n) = ( n*mu*pai(n) - (n+1)*pai(n-1) );
end

for n = 1:N
    
    % P = legendre(n, cos(theta)); Pn1 = P(2,1);  
    % P is an (n+1)x1 column vector
    % P(k,1) gives P_n^(k-1)[cos(theta)] for k in 1,2...(n+1)
    
    % dP = legendre_derivative(n, cos(theta)); dPn1 = dP(2,1);
    % Similar argument as above for derivative of assoc. leg. poly. 
    
    coeff = ((-1j)^n)*(2*n+1)/(n*(n+1));
    
    AA = A(n,x,m); BB = B(n,x,m); DD = Dh1(n,rho); HH = h1(n,rho); %Intermediate vars
    
    sumr = sumr + ( coeff*1j*AA*n*(n+1)*(HH/(rho))*cos(phi)*pai(n)*sin(theta) );
    
    sumt = sumt + ( coeff*cos(phi)*( (1j*AA*DD*tau(n)) + (BB*HH*pai(n)) ) );
    
    sump = sump - ( coeff*sin(phi)*( (1j*AA*DD*pai(n)) + (BB*HH*tau(n)) ) );
    
end

y = 4*pi*(r^2)*( abs(sumr)^2 + abs(sumt)^2 + abs(sump)^2 ); % RCS
%y = sumt;

% Er = sumr; Et = sumt; Ep = sump;
% Esph = [Er;Et;Ep];
% Ecart = sph2cartvec(Esph,90-(180*theta/pi),phi);
% Ex = abs(Ecart(1)); Ey = abs(Ecart(2)); Ez = abs(Ecart(3));


    function y = B(n,x,m)                       % scattered field coeff b_n
        y = ( psi(n,m*x)*Dpsi(n,x) - m*psi(n,x)*Dpsi(n,m*x) )...
            /( psi(n,m*x)*Dxi(n,x) - m*xi(n,x)*Dpsi(n,m*x) );
    end
        
    function y = A(n,x,m)                       % scattered field coeff a_n
        y = ( m*psi(n,m*x)*Dpsi(n,x) - psi(n,x)*Dpsi(n,m*x) )...
            /( m*psi(n,m*x)*Dxi(n,x) - xi(n,x)*Dpsi(n,m*x) );
    end

    function y = Dpsi(n,x)                     % Asymptotic derivative of riccati-bessel psi  
          %y = spj(n,x) + x*(n*spj(n-1,x) - (n+1)*spj(n+1,x) )/(2*n+1);
          y = ( spj(n,x) + x*(spj(n-1,x) - spj(n+1,x)) )/2;
    end

    function y = Dxi(n,x)                      % Asymptotic derivative of riccati-bessel xi
        
          %y = h1(n,x) + x*(n*h1(n-1,x) - (n+1)*h1(n+1,x) )/(2*n+1);
          y = ( h1(n,x) + x*(h1(n-1,x) - h1(n+1,x)) )/2;
    end

    function y = Dh1(n,x)                      % Asymptotic derivative (1/x)d(x h1(n,x))/dx
        
        %y =  (h1(n,x)/x) + (n*h1(n-1,x) - (n+1)*h1(n+1,x))/(2*n+1);
        y = ( (h1(n,x)/x) + h1(n-1,x) - h1(n+1,x) )/2;
    end

    function y = psi(n,x)                      % Riccati-Bessel psi(x) = xj(x)
        y = x*spj(n,x);
    end

    function y = xi(n,x)                       % Riccati-Bessel xi(x) = xh^(1)(x)
        y = x*h1(n,x);
    end

    function y = h1(n,x)                       % Spherical Bessel h^(1)_n(x)
        y = spj(n,x) - 1j*spy(n,x);   % Note the opposite sign convention (+j -> -j)
    end

    function y = spj(n,x)                      % Spherical Bessel j_n(x)
        y = sqrt(pi/(2*x))*besselj(n+0.5,x);
    end

    function y = spy(n,x)                      % Spherical Bessel y_n(x)
        y = sqrt(pi/(2*x))*bessely(n+0.5,x);
    end

end
