function [hsrp,GEV]=GPDfitwaves_funcion(hs_obs,time_obs,kind,rp);
  
  % kind = 1 = Annual Maxima - Generalised Extreme Value Distribution
  % kind = 2 = Peaks Over a Threshold - General Pareto Distribution 
  % rp   = return period (10, 50, 100 years, or else)

  %PF = it will Plot a Figure 
  %load wave_height_banks_peninsula.mat
  %
  hs_obs = hs_obs(:); time_obs = time_obs(:);

  if kind==1 % AM-GEV

    %cdfplot(hs_obs) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
    %return
    
    %tolerance =0.0001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
    %tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
    %ythr = GPthres(time_obs,hs_obs+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
    %disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m'])
    ARI = [1 2 5 10 20 50 100]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like. 
    %
    %[GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(hs_obs+tolerance,time_obs,tthr,ythr,ARI); % PF
    %
    %figure; clf
    %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
    %semilogx(ARI_POT,POT_sorted,'ok'), 
    %hold on
    %semilogx(ARI,GPD(1,:),'r',ARI,GPD(2,:),'--r',ARI,GPD(3,:),'--r')
    %title('GPD fitted to peaks over threshold'), xlabel('Average recurrence interval (years)'), ylabel('H_s')
    
    % Annual Maxima method
    %ARI = [1:100]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like. 
    AMy = AMcalc(time_obs,hs_obs,1); AM = AMy(:,3);
    [gev,gevci] = gevfit(AM);
    GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
    GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
    GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));

    [dif loc]=nanmin(abs(ARI-rp));
    hsrp=GEV(1,loc);
    
    %figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on
    %semilogx(ARI,GEV(1,:),'k','LineWidth',2)
    %semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
    %hold off
    %title('GEV fitted to annual maxima'), xlabel('Average recurrence interval  (years)'), ylabel('H_s')
    %
    %figure(2) 
    %hold on
    %semilogx(ARI,GEV(1,:),'b'), semilogx(ARI,GEV(2,:),'--b',ARI,GEV(3,:),'--b'), 
    %semilogx(AMari,AMsorted,'+b')
    %hold off
    %
    %legend('GPD-POT','GEV-AM');

  elseif kind==2 
    cdfplot(hs_obs) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
    %return
    
    tolerance =0.0001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
    tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
    ythr = GPthres(time_obs,hs_obs+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
    disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m'])
    ARI = [1 2 5 10 20 50 100]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like. 
    
    [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(hs_obs+tolerance,time_obs,tthr,ythr,ARI); % PF
    
    figure; clf
    semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
    semilogx(ARI_POT,POT_sorted,'ok'), 
    hold on
    semilogx(ARI,GPD(1,:),'r',ARI,GPD(2,:),'--r',ARI,GPD(3,:),'--r')
    title('GPD fitted to peaks over threshold'), xlabel('Average recurrence interval (years)'), ylabel('H_s')
    
    % Compare to annual maxima method. Don't need this but sometimes useful to 
    AMy = AMcalc(time_obs,hs_obs,1); AM = AMy(:,3);
    [gev,gevci] = gevfit(AM);
    GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
    GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
    GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
    
    figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on
    semilogx(ARI,GEV(1,:),'k','LineWidth',2)
    semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
    hold off
    title('GEV fitted to annual maxima'), xlabel('Average recurrence interval  (years)'), ylabel('H_s')
    
    figure(2) 
    hold on
    semilogx(ARI,GEV(1,:),'b'), semilogx(ARI,GEV(2,:),'--b',ARI,GEV(3,:),'--b'), 
    semilogx(AMari,AMsorted,'+b')
    hold off
    
    legend('GPD-POT','GEV-AM');

  end
