function y = RCS_huygen(dist_fac,thetap,phip) %[Ex,Ey,Ez]

load('tang.mat');
load('tang_phi.mat');
load('dErdphi.mat');
load('dErdth.mat');
%Variables: Ethet_1p48a, Ethet_1p5a, Ethet_1p52a
%Variables: Ephi_1p48a, Ephi_1p5a, Ephi_1p52a
%Variables: dErbydth, dErbydphi



thet = 0.0175:0.0175:3.1241;
phi = 0:0.0175:6.2657;

lambda=1; a=lambda; k = (2*pi)/lambda;

rp = dist_fac*a;

xp=rp*sin(thetap)*cos(phip);
yp=rp*sin(thetap)*sin(phip);
zp=rp*cos(thetap);

rp_vec = [xp;yp;zp];

vel = zeros(3,1);
tang = zeros(3,1); vell = zeros(3,1);

M = [(yp^2 + zp^2),-xp*yp,-xp*zp;...
      -xp*yp,(xp^2 + zp^2),-yp*zp;...
      -xp*zp,-yp*zp,(xp^2 + yp^2)];

for i = 1:length(thet)
    for j = 1:length(phi)

        xx = 1.5*a*sin(thet(i))*cos(phi(j)); yy = 1.5*a*sin(thet(i))*sin(phi(j));
        zz = 1.5*a*cos(thet(i));
        r_vec = [xx;yy;zz];

        dEphi_by_dr = (Ephi_1p52a(i,j) - Ephi_1p48a(i,j))/(0.00002*a);
        dEth_by_dr = (Ethet_1p52a(i,j) - Ethet_1p48a(i,j))/(0.00002*a);

        Ephi_by_r = Ephi_1p5a(i,j)/(1.5*a);
        Eth_by_r = Ethet_1p5a(i,j)/(1.5*a);
        Ep = 1.5*a*Ephi_by_r;
        Et = 1.5*a*Eth_by_r;

        dEr_by_dt = dErbydth(i,j)/(1.5*a);
        dEr_by_dp = dErbydphi(i,j)/(1.5*a);

        sint = sin(thet(i)); cost = cos(thet(i));
        sinp = sin(phi(j)); cosp = cos(phi(j));

        expo = exp(-1j*k*norm(r_vec - rp_vec));

        vel(1) = vel(1) + ( sint*sinp*(dEphi_by_dr + Ephi_by_r - dEr_by_dp/sint) ...
                           -cosp*cost*sint*(dEth_by_dr + Eth_by_r - dEr_by_dt) )*expo;

        vel(2) = vel(2) - ( sint*cosp*(dEphi_by_dr + Ephi_by_r - dEr_by_dp/sint) ...
                            +sinp*cost*sint*(dEth_by_dr + Eth_by_r - dEr_by_dt) )*expo;

        vel(3) = vel(3) + sint^2*(dEth_by_dr + Eth_by_r - dEr_by_dt)*expo;


        tang(1) = -( Et*sinp + Ep*cost*cosp );
        tang(2) =  ( Et*cosp - Ep*cost*sinp );
        tang(3) =  Ep*sint;

        vell(1) = vell(1) + ((yp - yy)*tang(3) - (zp - zz)*tang(2))*sint*expo;
        vell(2) = vell(2) + ((zp - zz)*tang(1) - (xp - xx)*tang(3))*sint*expo;
        vell(3) = vell(3) + ((xp - xx)*tang(2) - (yp - yy)*tang(1))*sint*expo;
    end
end

terma = (0.0175^2 * (1.5*a)^2 / (4*pi*rp^3)) * M*vel;

termb = (-1j*k*0.0175^2 *(1.5*a)^2 / (4*pi*rp^2)) * vell;

Efar = terma + termb;

y = 4*pi*rp^2 * ( abs(Efar(1))^2 + abs(Efar(2))^2 + abs(Efar(3))^2 );
% Ex = abs(Efar(1));
% Ey = abs(Efar(2));
% Ez = abs(Efar(3));

end
