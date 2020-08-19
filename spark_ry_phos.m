colors = repmat('krgbmc',1,500) ;

dt = 1e-5 ;        % ms, adjust so that transition probability <= 0.01
dt_record = 0.1 ;

% % Protocol: open single RyR at 'interval'
% % Then continue to run for 'timeafter' milliseconds
interval = 0 ;
%timeafter = 500 ;
t_end = 1000; 

trials = 1 ;    % Number of release units to simulate

F = 96.485 ;       % Faraday's constant, C/mmol
% %% Needed to convert from flux to current

V_ds = 1.0000e-12 ;
V_JSR = 1.6000e-12 ;

tau_efflux = 1.78e-3 ;           % ms
% % Per Ramay paper, DCa = 225 (um)^2/s, dx = 0.02 um
tau_refill = 6.55 ;
% % Fits total JSR [Ca] recovering with tau = 90 ms

EJequiv = 0.1 ;
% % Defines coupling between RyRs

D_ryr = 2.2e-12 ;                    % ms-1
% % Defines current through single RyR = 0.42 pA

% %% RyR gating parameters
kr_minus = 0.48 ;              % ms-1 (mean open time = 2 ms) 
kr_plus_max = 42 ;             % ms-1

kr_plus_max_ryanodine = 23052865.2;
Km_r_max_orig = 28.17 ;             % uM
alpha_r_orig = 9.3e-3 ;             % unitless
hill = 4 ;                     % exponent
CaJSR_diastolic = 1000;
Km_fixed = Km_r_max_orig - alpha_r_orig*CaJSR_diastolic;

alpha_r = 1.0e-3;
Km_r_max = Km_fixed + alpha_r*CaJSR_diastolic;

N_RyR_normal = 27;
N_RyR_ryanodine = 1;
N_RyR = N_RyR_normal + N_RyR_ryanodine; 
kcoup = exp(2*EJequiv/(N_RyR-1)) ;

k_PKA_RyR = 1.9; %TODO: need to test for steady state value
PKAC_I = 1; %TODO: this is not a known value
PPI = 0.89; %TODO: this is from the saucerman paper
Km_PKA_RyR = 50; %TODO: is this a known value?
k_PPI_RyR = 0.3 ; %TODO: is this a known value?
Km_PPI_RyR = 35; %TODO: is this a known value?

% %% Buffering parameters
% Subspace buffers, calmodulin, SR sites, SL sites
bt =  [24 47 900] ;            % uM
n_buffers = length(bt) ;
% Factor of 1e-3 to convert from s-1 to ms-1
kp =  1e-3*[100 115 115] ;     % uM^-1 ms^-1
km =  1e-3*[38 100 1000] ;     % ms^-1
% JSR buffer CSQ
CSQ = 30e3 ;       % uM
KCSQ = 630 ;       % uM

% %% Fixed ionic concentrations
Camyo = 0.1 ;
CaNSR = 1000 ;

% %% Initialize arrays to hold results
filename = 'data.mat' ;

%t_end = interval + timeafter ;
%iterations = round(t_end/dt) ;
outputs = round(t_end/dt) ;
plottime = 0:dt_record:(outputs-1)*dt_record ;


Cads_all = zeros(outputs,trials) ;
CaJSR_all = zeros(outputs,trials) ;
%CaJSR_all = [];
Irel_all = zeros(outputs,trials) ;
Nopen_all = zeros(outputs,trials) ;
%Nopen_all = [];
Nopen_phos_all = zeros(outputs, trials);

tic ;
  
for j=1:trials

    %% Initial conditions
    Cads = Camyo ;
    CaJSR = CaNSR ;
    nopen_phos = 0 ;
    nopen_dephos = 0; %ADDITION
    nopen_ryanodine = 0;
    nopen_normal = 0;
    nopen = nopen_phos + nopen_dephos;
    b = bt.*(km./kp)./(km./kp+Cads);
    writedex = 1 ;

    
    Km_r_ryanodine = Km_r_max - alpha_r*CaJSR;
    
     neverspark = 1 ;
    
    tlast = -dt ;
    time = 0;
    
    %i = 1;

    while time < t_end


      time = tlast + dt; 
      if (time >= interval && tlast < interval)
        nopen = nopen + 5 ;
      end
      if (time >= interval+10.0 && nopen < 1 && neverspark)
        break ;
      end
      
      nopen = nopen_phos + nopen_dephos;
      nclosed_ryanodine = N_RyR_ryanodine - nopen_ryanodine;
      nclosed_normal = (N_RyR - N_RyR_ryanodine) - nopen_normal;
      
        %plb_p shoulld be the # of n_phos
      dRyRp_dt = - ((k_PPI_RyR * PPI * N_RyR_phos) / (Km_PPI_RyR + N_RyR_phos)) + ((k_PKA_RyR * PKAC_I * N_RyR_dephos) / (Km_PKA_RyR + N_RyR_dephos));
      if (N_RyR_phos + dRyRp_dt < N_RyR)
         N_RyR_phos = N_RyR_phos + dRyRp_dt;
      end
      N_RyR_dephos = N_RyR - N_RyR_phos;
      
      J_ryr = nopen*D_ryr*(CaJSR-Cads)/V_ds ;    % uM/ms
      I_ryr = 1e6*J_ryr*2*F*V_ds ;               % pA
      J_efflux = (Camyo-Cads)/tau_efflux ;
      J_refill = (CaNSR-CaJSR)/tau_refill ;

      db_dt = -kp.*b*Cads + km.*(bt-b) ;
      
      J_buff = sum(db_dt) ;
      
      B_JSR = (1 + CSQ*KCSQ/(KCSQ+CaJSR)^2)^-1 ;
      
   
      if (mod(time,dt_record) == 0)
        %disp(time)
        Cads_all(writedex,iiii) = mean(Cads) ;
        CaJSR_all(writedex,iiii) = mean(CaJSR) ;
        Irel_all(writedex,iiii) = I_ryr ;
        Nopen_phos_all(writedex,iiii) = nopen_phos ;
        Nopen_dephos_all(writedex,iiii) = nopen_dephos;
        Nopen_all(writedex, iiii) = nopen_phos + nopen_dephos;
        writedex = writedex + 1 ;
      end
%       
      
      Km_r = Km_r_max - alpha_r*CaJSR ;
      kr_plus = kr_plus_max*Cads^hill/(Cads^hill + Km_r^hill) ;
      kr_plus_phos = kr_plus_max_phos *Cads^hill/(Cads^hill + Km_r^hill) ;
      kr_plus_ryanodine = kr_plus_max_ryanodine * 0.1 ^ hill / 0.1 ^ hill + Km_r_ryanodine ^ hill;
      
      % ADDITION
      pincrease_phos = dt*(N_RyR_phos - nopen_phos)*kr_plus_phos*kcoup^(2*nopen_phos + 1 - N_RyR) ;
      pincrease_dephos = dt*(N_RyR_dephos - nopen_dephos)*kr_plus*kcoup^(2*nopen_dephos + 1 - N_RyR) ;
      pdecrease_phos = dt*nopen_phos*kr_minus_phos*kcoup^(2*(N_RyR_phos - nopen_phos) + 1 - N_RyR) ;
      pdecrease_dephos = dt*nopen_dephos*kr_minus*kcoup^(2*(N_RyR_dephos - nopen_dephos) + 1 - N_RyR) ;
      pincrease_ryanodine = dt*nclosed_ryanodine*kr_plus_ryanodine*kcoup^(2*nopen + 1 - N_RyR) ;
      
      if (rand < pincrease_phos)
          if (nopen_phos < N_RyR_phos - 1 && N_RyR_phos < N_RyR);
            nopen_phos = nopen_phos + 1 ;
          end
      end
      
      if (rand < pdecrease_phos)
        nopen_phos = nopen_phos - 1 ;
      end
      
      if (rand < pincrease_dephos)
         if (nopen_dephos < N_RyR_dephos - 1 && N_RyR_dephos < N_RyR)
            nopen_dephos = nopen_dephos + 1;
         end
      end
      if (rand < pdecrease_dephos)
          nopen_dephos = nopen_dephos - 1;
        
      end

      if ((nopen) >= 5) 
        neverspark = 0 ;
      end
      
      
      
      if(rand < pincrease_ryanodine)
          if(nopen + 1 <= N_RyR && N_RyR_ryanodine + 1 <= N_RyR)
              nopen_ryanodine = nopen_ryanodine + 1;
              nopen_phos = nopen_phos + 1;
              N_RyR_ryanodine = N_RyR_ryanodine + 1;
              N_RyR_normal = N_RyR - N_RyR_ryanodine;
          end 
      end
      
      nopen = nopen + nopen_ryanodine;
      nopen_normal = nopen - nopen_ryanodine;
      % nopen should also be equal to nopen_phos + nopen_dephos
      
      dCads_dt = J_efflux + J_d + J_ryr + J_buff ;
      dCaJSR_dt = B_JSR *(J_refill - J_ryr*V_ds/V_JSR) ;     

      Cads = Cads + dt*dCads_dt ;
      CaJSR = CaJSR + dt*dCaJSR_dt ;
      
      b = b + dt*db_dt ;
      
      tlast = time ;
%       plottime(i) = time; 
%       
%       Cads_all(i,j) = Cads;
%       CaJSR_all(i,j) = CaJSR ;
%       Irel_all(i,j) = 1e6*J_ryr*2*F*V_ds ;
%       Nopen_all(i,j) = nopen ;
%       
%       i = i + 1;
    end
    
     %% write values at last time point t_end after loop finished
    Cads_all(end,iiii) = Cads ;
    CaJSR_all(end,iiii) = CaJSR ;
    
    Irel_all(end,iiii) = 1e6*J_ryr*2*F*V_ds ;
    Nopen_all(end,iiii) = nopen ;
    Nopen_phos_all(end,iiii) = nopen_phos;
   

end

%save(filename,'plottime','Cads_all','CaJSR_all','Irel_all','Nopen_all') ;
    
minutes = toc/60 ;
disp(['Simulation completed in ',num2str(minutes),' minutes'])

figure 
hold on
plot(plottime, Nopen_all, 'k');
plot(plottime, Nopen_phos_all, 'b');
xlabel('Time in ms');
ylabel('Number of RyRs');
legend('Open receptors', 'Phosphorylated');
hold off

figure 
plot(plottime, Cads_all);
xlabel('Time in ms');
ylabel('[Ca2+] in subspace');
% sparkdices = find(max(Nopen_all) > 5) ;
%
% figure
% plot(plottime,CaJSR_all(:,sparkdices))
% hold on 
% 
% CaJSR_avg = mean(CaJSR_all(:,sparkdices),2) ;
% plot(plottime,CaJSR_avg,'b','LineWidth',2.25)
% print -depsc freeCaJSR
% 
% CaJSRtot_all = CaJSR_avg + ...
%   CaJSR_avg*CSQ./(KCSQ + CaJSR_avg) ;
% CaJSRtot_all = CaJSRtot_all/max(max(CaJSRtot_all)) ;
% figure
% plot(plottime,CaJSRtot_all)
% hold on
% plot(plottime, 1 - exp(-plottime/90),'r','LineWidth',2.5)
% print -depsc totalCaJSR

