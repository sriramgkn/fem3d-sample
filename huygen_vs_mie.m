dist_fac=100; phip=0; N=15; lambda=1; epsi=2.56;

a = lambda;
% 
% thetp = 0:0.04:pi; % theta prime, at far field
% 
% y = zeros(length(thetp),1);
% yy = zeros(length(thetp),1); %For huygen
% yyy = zeros(1,length(thetp));
% 
% rcs_db = zeros(1,length(thetp));
% rcs_db_2 = zeros(1,length(thetp)); %For huygen
% rcs_db_3 = zeros(1,length(thetp));
% 
% for i=1:length(thetp)
%     y(i) = RCS(dist_fac,thetp(i),phip,N,lambda,epsi,a);
%     %yyy(i) = RCS(1.5,thetp(i),phip,N,lambda,epsi,a);
%     
%     %epsi=2.56, lambda=1, a=lambda (implicit)
%     yy(i) = RCS_huygen(dist_fac,thetp(i),phip); 
%     
%     
%     rcs_db(i) = 10 * log10( y(i)/(lambda^2) );
%     rcs_db_2(i) = 10 * log10( yy(i)/(lambda^2) );
%     %rcs_db_3(i) = 10 * log10( yyy(i)/(lambda^2) );
% end
% 
% plot(180*thetp/pi,rcs_db,180*thetp/pi,rcs_db_2,'LineWidth',1.8);

load('tang.mat');
load('tang_phi.mat');
load('dErdphi.mat');
load('dErdth.mat');
thet = 0.0175:0.0175:3.1241;
phi = 0:0.0175:6.2657;

for i = 1:length(thet)
    dEphi_by_dr = (Ephi_1p52a(i,1) - Ephi_1p48a(i,1))/(0.00002*a);
    dEth_by_dr = (Ethet_1p52a(i,1) - Ethet_1p48a(i,1))/(0.00002*a);
    Ephi_by_r = Ephi_1p5a(i,1)/(1.5*a);
    Eth_by_r = Ethet_1p5a(i,1)/(1.5*a);
    dEr_by_dt = dErbydth(i,1)/(1.5*a);
    dEr_by_dp = dErbydphi(i,1)/(1.5*a);
    sint = sin(thet(i)); cost = cos(thet(i));
    sinp = sin(phi(1)); cosp = cos(phi(1));
    
    term1 = dEth_by_dr + Eth_by_r - dEr_by_dt;
    term2 = dEphi_by_dr + Ephi_by_r - dEr_by_dp/sint;
    % (r,th,phi): (0, -term1, -term2)
    v_pol = [0;-term1;-term2]; % r x (nabla x E) in sph polar
    U = [sint*cosp, sint*sinp, cost;...
         cost*cosp, cost*sinp, -sint;...
         -sinp, cosp, 0];
     
    v_cart = U'*v_pol; % r x (nabla x E) in cartesian
    term2
end

% subplot(1,3,1);
% plot(180*thetp/pi,y(:,1),180*thetp/pi,yy(:,1));
% legend('Mie Series','Huygen');
% ylabel('E_x'); xlabel('\theta');
% 
% subplot(1,3,2);
% plot(180*thetp/pi,y(:,2),180*thetp/pi,yy(:,2));
% legend('Mie Series','Huygen');
% ylabel('E_y'); xlabel('\theta');
% 
% subplot(1,3,3);
% plot(180*thetp/pi,y(:,3),180*thetp/pi,yy(:,3));
% legend('Mie Series','Huygen');
% ylabel('E_z'); xlabel('\theta');
% 
% suptitle('X, Y, Z components of far-field vs \theta');

legend('Mie Series','Huygen');
ylabel('10 log(\sigma/\lambda^2)'); 
xlabel('Elevation angle in degrees (\theta)'); 
title('\sigma (dB) vs \theta: \epsilon_r= 2.56, \lambda=1m, r=100a, \phi=0, a = \lambda');
