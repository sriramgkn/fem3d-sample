
N=15; lambda=1; epsi=2.56;

% Array=csvread('book_mie.csv');
% col1 = Array(:, 1);
% col2 = Array(:, 2);

a = lambda;
% a2 = 0.5*lambda;
% a3 = 0.2*lambda;

thet = 0.0175:0.0175:3.1241;
phi = 0:0.0175:6.2657;

%Variables: Ethet_1p48a, Ethet_1p5a, Ethet_1p52a
%Variables: Ephi_1p48a, Ephi_1p5a, Ephi_1p52a

% Ethet_1p48a = zeros(length(thet),length(phi));
% Ethet_1p5a = zeros(length(thet),length(phi));
% Ethet_1p52a = zeros(length(thet),length(phi));

% Ephi_1p48a = zeros(length(thet),length(phi));
% Ephi_1p5a = zeros(length(thet),length(phi));
% Ephi_1p52a = zeros(length(thet),length(phi));

%dErbydphi = zeros(length(thet),length(phi));

dErbydth = zeros(length(thet),length(phi));

% y2 = zeros(1,length(thet));
% y3 = zeros(1,length(thet));
% rcs_db = zeros(1,length(thet));
% rcs_db_2 = zeros(1,length(thet));
% rcs_db_3 = zeros(1,length(thet));

temp2 = RCS(1.5,0,0,N,lambda,epsi,a);

for j=1:length(phi)
    
%     temp = RCS(1.5,thet(i),0,N,lambda,epsi,a);
      temp = temp2;
    
    for i=1:length(thet)
        
%         if j==length(phi)
%             dErbydphi(i,j) = (RCS(1.5,thet(i),6.2831,N,lambda,epsi,a)-temp)/0.0175;
%         else
%             dErbydphi(i,j) = (RCS(1.5,thet(i),phi(j+1),N,lambda,epsi,a)-temp)/0.0175;
%             temp = temp + 0.0175*dErbydphi(i,j);
%         end
            

          dErbydth(i,j) = (RCS(1.5,thet(i),phi(j),N,lambda,epsi,a)-temp)/0.0175;
          temp = temp + 0.0175*dErbydth(i,j);


        
%           Ethet_1p48a(i,j) = RCS(1.49999,thet(i),phi(j),N,lambda,epsi,a);
%           Ethet_1p5a(i,j) = RCS(1.5,thet(i),phi(j),N,lambda,epsi,a);
%           Ethet_1p52a(i,j) = RCS(1.50001,thet(i),phi(j),N,lambda,epsi,a);
%           
%           Ephi_1p48a(i,j) = RCS(1.49999,thet(i),phi(j),N,lambda,epsi,a);
%           Ephi_1p5a(i,j) = RCS(1.5,thet(i),phi(j),N,lambda,epsi,a);
%           Ephi_1p52a(i,j) = RCS(1.50001,thet(i),phi(j),N,lambda,epsi,a);
          
          % y2(i) = RCS(dist_fac,thet(i),phi,N,lambda,epsi,a2);
          % y3(i) = RCS(dist_fac,thet(i),phi,N,lambda,epsi,a3);
          
          % rcs_db(i) = 10 * log10( y(i)/(lambda^2) );
          % rcs_db_2(i) = 10 * log10( y2(i)/(lambda^2) );
          % rcs_db_3(i) = 10 * log10( y3(i)/(lambda^2) );
    end
end

% save('tang.mat','Ethet_1p48a','Ethet_1p5a','Ethet_1p52a');

save('dErdth.mat','dErbydth');

% plot(180*thet/pi,y); 
% ylabel('4\pi r^2|E_{scat}|^2'); 
% xlabel('Elevation angle in degrees (\theta)'); 
% title('RCS vs \theta: \epsilon_r=10.0, \lambda=30cm, a=0.2*\lambda,
% r=10a, \phi=\pi');

% plot(180*thet/pi,rcs_db,180*thet/pi,rcs_db_2,180*thet/pi,rcs_db_3,'linewidth',1.2); 
% 
% legend('a = \lambda','a = 0.5\lambda','a = 0.2\lambda');
% ylabel('10 log(\sigma/\lambda^2)'); 
% xlabel('Elevation angle in degrees (\theta)'); 
% title('\sigma (dB) vs \theta: \epsilon_r= 2.56, \lambda=1m, r=100a, \phi=0');
% 
% figure; plot(180*thet/pi,rcs_db,col1,col2,'linewidth',1.2);
% legend('Our code','Textbook');
% ylabel('10 log(\sigma/\lambda^2)'); 
% xlabel('Elevation angle in degrees (\theta)'); 
% title('\sigma (dB) vs \theta: \epsilon_r= 2.56, \lambda=1m, a=\lambda, r=100a, \phi=0');


