function [corr_data]=bias_correction(obsi,modi,orig_data,plot_data)
  
  % obsi = reference data (satellite or in situ obs) 
  % modi = model interpolated and cleaned to match obsi
  % orig_data = model original data that will be used in correction
  % corr_data = final bias corrected data
  % plot_data = 0 no plot; 1 =plot 

  %% Reference
  %[1] Singh, A.; Gaurav, K.; Meena, G.K.; Kumar, S. Estimation of Soil Moisture Applying Modified Dubois Model to Sentinel-1; A Regional Study from Central India. Remote Sens. 2020, 12, 2266.
  %[2] Reichle, R.H.; Koster, R.D. Bias reduction in short records of satellite soil moisture. Geophys. Res. Lett. 2004, 31, L19501.


  Ref_data=obsi(:);
  Biased_data=modi(:);

  bias=sum(Biased_data-Ref_data)/length(Biased_data); %Bias value before bais correction
  %% CDF Matching (Eq. 1) of [2]
  Prob_ref=((1:length(Ref_data))./length(Ref_data))';
  Prob_biased = Prob_ref;
  Prob_corr=((1:length(orig_data))./length(orig_data))';
  Biased_interpolated=sort(Biased_data); 

  n=3; %degree of the polynomial
  p= polyfit(Biased_interpolated,sort(Ref_data)-Biased_interpolated, n); % Polynomial fitting between the biased data and difference of referenced and interpolated data.

  Corrections= polyval(p,Biased_data); % Evaluates the polynomial p at each point in Biased_data plus the Biased_data
  Corrected_Data= Corrections + Biased_data; % Evaluates the polynomial p at each point in Biased_data plus the Biased_data
  
  curve_fitted = polyval(p,Biased_interpolated); % Evaluates the polynomial p at each point in Biased_data plus the Biased_data

  corr_data = polyval(p,orig_data) + orig_data; % Evaluates the polynomial p at each the original timeseires

  bias_after=sum(Corrected_Data-Ref_data)/length(Corrected_Data); %Bias value after bais correction

  if plot_data==1

    %% Plotting
    figure
    plot(sort(Ref_data), Prob_ref , 'k-', 'linewidth',1.5)
    hold on
    plot(sort(Biased_data),Prob_biased, 'k--','linewidth',1.5)
    hold on
    plot( sort(Corrected_Data),Prob_data(Corrected_Data), 'k:', 'linewidth',1.5)
    xlabel('Variable')
    ylabel('CDF')
    leg=legend ('Reference data','Biased data','Corrected data');
    set(leg,'location','best')
    set(gcf,'color','white')

    figure; hold on
    plot(Biased_interpolated,sort(Ref_data)-Biased_interpolated,'b')
    plot(Biased_interpolated,curve_fitted,'k')
    xlabel('Variable')
    ylabel('Corrections')
    leg=legend ('Reference - Model','Corrections');
    set(leg,'location','best')
    set(gcf,'color','white')

  end

end

