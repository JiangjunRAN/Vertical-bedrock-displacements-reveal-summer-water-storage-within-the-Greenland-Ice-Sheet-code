%% estimate meltwater storage from SMB and GNSS
function [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig,formal_error_calibrated_E_r]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r(npol,year_reg,invbeta_reg,ErMeanReg,nit,damping_scaling,time,GNSS_loading_time_series,SMB_time_series,liquid_p_time_series,runoff_time_series,year_beg,year_end)
%***********************************************************************************************%
% Copyright (c) 2024-2099 by Gravity and Cryosphere Lab @ SUSTecch, all rights reserved.
% Author: Jiangjun Ran
% E-mail: ranjj@sustech.edu.cn
% cite: Ran, J., Ditmar, P., van den Broeke, M. R., et al. (in press). Vertical bedrock shifts reveal summer water storage in Greenland ice sheet. Nature.
%***********************************************************************************************%
%
% Input parameters:
% npol: Number of terms in the polynomial describing the ice discharge (e.g., 2 for a linear function)
% year_reg: (when > 0) the year in which 'E_r' is forced to be nearly -1 (i.e., the water mass is forced to be nearly zero)
% invbeta_reg: (when > 0) the predefined value of parameter '1/beta'
% ErMeanReg: (when > 0) the mean E_r is forced to be zero
% nit: number of iterations
% damping_scaling: scaling to be applied to the damping factor at each iteration (1. if no damping is to be applied)
% ...
%
    % initialization
    nia=npol;
%     year_beg=2009;
%     year_end=2015;
    ny=year_end-year_beg+1
    n_unknowns=nia+ny+2
    k_year_reg=year_reg-year_beg+1
%   x_0(1:5,1)=0.001;
    x_0(1:nia,1)=0.0;
    x_0(nia+1:nia+ny,1)=0.0;
    x_0(nia+ny+1,1)=1.0;
    if(year_reg>=year_beg&&year_reg<=year_end)
      x_0(nia+k_year_reg,1)=-1.0; % Setting the initial meltwater mass to be 0 in the regularized year
    end
% Setting the initial 1/beta
    x_0(nia+ny+2,1)=0.1; 
    if invbeta_reg>0
      x_0(nia+ny+2,1)=invbeta_reg;
    end
%   
    delta_x(1:n_unknowns,1)=0.0;
    DesignMatrix(1:length(time),1:n_unknowns)=0.0;
    DesignMatrix2(1:length(time),1:n_unknowns)=0.0;
    intergated_upstream_runoff_per_year_orig(1:length(time),1:ny)=0.0;
    delta_tao=1.0/365.25;
    damping=1.;
    n_pseudo=0;
    if(year_reg>=year_beg&&year_reg<=year_end)
      n_pseudo=n_pseudo+1;
    end
    if invbeta_reg>0
      n_pseudo=n_pseudo+1;
    end
    if ErMeanReg>0
      n_pseudo=n_pseudo+1;
    end
    if n_pseudo>0
      DesignMatrix_pseudo(1:n_pseudo,1:n_unknowns)=0.0;
      weight_pseudo=1.e10;
  % design matrix of pseudo-observations (input for regularization)
      k_pseudo=0;
      if(year_reg>=year_beg&&year_reg<=year_end)
        k_pseudo=k_pseudo+1;
        DesignMatrix_pseudo(k_pseudo,nia+k_year_reg)=1.0;
      end
      if invbeta_reg>0
        k_pseudo=k_pseudo+1;
        DesignMatrix_pseudo(k_pseudo,nia+ny+2)=1.0;
      end
      if ErMeanReg>0
        k_pseudo=k_pseudo+1;
        DesignMatrix_pseudo(k_pseudo,nia+1:nia+ny)=1.0/ny;
      end
  % regularization matrix
      RegMatrix=weight_pseudo*DesignMatrix_pseudo'*DesignMatrix_pseudo;
    end
  % design matrix
    time_rel=time-mean(time);
    counter=1;
    while (counter<=nit)
      DesignMatrix(1:length(time),1:n_unknowns)=0.0;
      DesignMatrix2(1:length(time),1:n_unknowns)=0.0;
        for k=1:length(time_rel)
            % Inter-annual variability
	    time_power=1.;
 	    for kpol=1:npol
              DesignMatrix(k,kpol)=time_power;
              DesignMatrix2(k,kpol)=time_power;
	      time_power=time_power*time_rel(k);
	    end
%            % epsilon_p-epsilon_r
%            temp_sum=0.0;
%            temp_sum2=0.0;
%            for n=1:k
%                temp_sum=temp_sum+liquid_p_time_series(n,1)*delta_tao; % ? liquid_p_time_series(n,1) has already been accumulated from the first epoch. It should be changed to daily water production?
%                temp_sum2=temp_sum2+(-1.0*runoff_time_series(n,1)*delta_tao); % daily runoff?
%            end
%            DesignMatrix(k,3)=temp_sum+temp_sum2;
            % epsilon_r, theta
            temp_sum(1:ny)=0.0;
            temp_sum2(1:ny)=0.0;
            for n=1:k
	      year=floor(time(n));
	      ky=year-year_beg+1;
	      if ky>=1 && ky<=ny
	        for ly=1:ny
		  runoff_cur=0.0;
		  if ly==ky
		    runoff_cur=runoff_time_series(n,1);
		  end
%                 temp_sum(ly)=temp_sum(ly)+runoff_cur*(exp(-1.0*(x_0(nia+ny+1,1))*(time_rel(k,1)-time_rel(n,1)))-1.0)*delta_tao; % daily runoff?
                  temp_sum(ly)=temp_sum(ly)+runoff_cur*(-1.0)*delta_tao; % daily runoff?
                  temp_sum2(ly)=temp_sum2(ly)+runoff_cur*exp(-1.0/(x_0(nia+ny+2,1))*(time_rel(k,1)-time_rel(n,1)))*delta_tao; % daily runoff?
		end
	      end
            end
            DesignMatrix(k,nia+ny+1)=0.0;
            DesignMatrix2(k,nia+ny+1)=0.0;
            for ly=1:ny
              DesignMatrix(k,nia+ly)=temp_sum(ly)+x_0(nia+ny+1,1)*temp_sum2(ly);
              DesignMatrix2(k,nia+ly)=DesignMatrix(k,nia+ly);
              DesignMatrix(k,nia+ny+1)=DesignMatrix(k,nia+ny+1)+(1.0+x_0(nia+ly,1))*temp_sum2(ly);
              DesignMatrix2(k,nia+ny+1)=DesignMatrix2(k,nia+ny+1)+temp_sum2(ly);
	      intergated_upstream_runoff_per_year_orig(k,ly)=temp_sum(ly);
	    end
%            % gamma
%            temp_sum2(1:ny)=0.0;
%            for n=1:k
%	      year=floor(time(n));
%	      ky=year-year_beg+1;
%	      if ky>=1 && ky<=ny
%	        temp_sum2(ky)=temp_sum2(ky)+runoff_time_series(n,1)*exp(-1.0/(x_0(nia+2*ny+ky,1))*(time_rel(k,1)-time_rel(n,1)))*delta_tao;
%	      end
%            end
%            for ky=1:ny
%              DesignMatrix(k,nia+ny+ky)=temp_sum2(ky);
%	    end
            % 1/beta
            temp_sum2(1:ny)=0.0;
            for n=1:k
	      year=floor(time(n));
	      ky=year-year_beg+1;
	      if ky>=1 && ky<=ny
 	        for ly=1:ny
		  runoff_cur=0.0;
		  if ly==ky
		    runoff_cur=runoff_time_series(n,1);
		  end
%                 temp_sum(ky)=temp_sum(ky)+(x_0(nia+ny+ky,1))*runoff_time_series(n,1)*exp(-1.0/(x_0(nia+2*ny+ky,1))*(time_rel(k,1)-time_rel(n,1)))*(time_rel(n,1)-time_rel(k,1))*delta_tao;
                  temp_sum2(ly)=temp_sum2(ly)+x_0(nia+ny+1,1)*(1.0+x_0(nia+ly,1))*runoff_cur*exp(-1.0/(x_0(nia+ny+2,1))*(time_rel(k,1)-time_rel(n,1)))*(time_rel(k,1)-time_rel(n,1))/x_0(nia+ny+2,1)/x_0(nia+ny+2,1)*delta_tao;
		end
	      end
            end
            DesignMatrix(k,nia+ny+2)=0.0;
            for ky=1:ny
              DesignMatrix(k,nia+ny+2)=DesignMatrix(k,nia+ny+2)+temp_sum2(ky);
	    end
        end
        nm=DesignMatrix'*DesignMatrix; % Normal matrix
	if (n_pseudo>0)
	  nm=nm+RegMatrix;		 % Regularization
	end
        cm=inv(nm);  % Covariance matrix
        time_series_to_fit=GNSS_loading_time_series-SMB_time_series-DesignMatrix2(:,1:nia+ny+1)*x_0(1:nia+ny+1,1);
        rhs=DesignMatrix'*time_series_to_fit; % Right-hand side vectors
        sol=cm*rhs;  % Solutions
        delta_x=sol*damping;
	damping=damping*damping_scaling
        x_new=x_0+delta_x;
        residual_postfit=time_series_to_fit-DesignMatrix*delta_x;
        sigma_postfit=sqrt((residual_postfit'*residual_postfit)/(length(residual_postfit)-length(delta_x)));
%	disp('CP-1');
%        disp(x_0);
	x_0=x_new;
%        for ky=1:ny
%          if x_0(nia+ky)<-0.9
%	    x_0(nia+ky)=-0.9;
%	  end
%        end
        if x_0(nia+ny+1)<0.01
	  x_0(nia+ny+1)=0.01;
	end
        if x_0(nia+ny+2)>2.
	  x_0(nia+ny+2)=2.;
	end
        if x_0(nia+ny+2)<0.001
	  x_0(nia+ny+2)=0.001;
	end
        disp(counter);
        disp([delta_x(1:nia,1) x_new(1:nia,1) x_0(1:nia,1)]);
        disp([delta_x(nia+1:nia+ny,1) x_new(nia+1:nia+ny,1) x_0(nia+1:nia+ny,1)]);
        disp([delta_x(nia+ny+1,1) x_new(nia+ny+1,1) x_0(nia+ny+1,1)]);
        disp([delta_x(nia+ny+2,1) x_new(nia+ny+2,1) x_0(nia+ny+2,1)]);
        disp([52.*delta_x(nia+ny+2,1) 52.*x_new(nia+ny+2,1) 52.*x_0(nia+ny+2,1)]);
%       disp([delta_x(nia+2*ny+1:nia+3*ny,1) x_new(nia+2*ny+1:nia+3*ny,1) x_0(nia+2*ny+1:nia+3*ny,1)]);
       counter=counter+1;
    end
    
    a=x_0(1,1);
    trend=x_0(2,1);
    c2=x_0(3,1);
    c3=x_0(4,1);
    E_r(1:ny,1)=x_0(nia+1:nia+ny,1);
 %   gamma(1:ny,1)=x_0(nia+ny+1:nia+2*ny,1);
    theta=x_0(nia+ny+1,1);
    invbeta=x_0(nia+ny+2,1);
    formal_error_calibrated_invbeta=sqrt(sigma_postfit^2*cm(nia+ny+2,nia+ny+2));
    % formal error E_r
    for n_meltwater=1:ny
        formal_error_calibrated_E_r(n_meltwater,1)=sqrt(sigma_postfit^2*cm(nia+n_meltwater,nia+n_meltwater));
    end
    time_series_simulated=DesignMatrix2(:,1:nia+ny+1)*x_0(1:nia+ny+1,1);
    discharge=DesignMatrix(:,1:npol)*x_0(1:npol,1);
    meltwater=DesignMatrix(:,nia+ny+1)*x_0(nia+ny+1,1);
    formal_error_calibrated_x_0_nia_ny_1=sqrt(sigma_postfit^2*cm(nia+ny+1,nia+ny+1));
    covariance_matrix_meltwarer=sqrt(DesignMatrix(:,nia+ny+1)*formal_error_calibrated_x_0_nia_ny_1^2*DesignMatrix(:,nia+ny+1)');
    for n_meltwater=1:length(DesignMatrix(:,nia+ny+1))
    calibrated_formal_error_meltwarer(n_meltwater,1)=covariance_matrix_meltwarer(n_meltwater,n_meltwater);
    end
    meltwater_orig(1:length(time),1)=DesignMatrix2(:,nia+ny+1)*x_0(nia+ny+1,1);
    for ly=1:ny
      intergated_upstream_runoff_per_year(1:length(time),ly)=intergated_upstream_runoff_per_year_orig(:,ly).*(1.+x_0(nia+ly,1));
    end

%    DM_col=DesignMatrix(:,nia+2);

%    eta(1:ny,1)=1.+E_r(1:ny,1);
%    k(1:ny,1)=eta(1:ny,1).*beta(1:ny,1)./gamma(1:ny,1);
%    Dt(1:ny,1)=(gamma(1:ny,1)-eta(1:ny,1))./(eta(1:ny,1).*beta(1:ny,1));
%    disp([k Dt]);

end
