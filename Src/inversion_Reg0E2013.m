%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
% Estimating the summer water storage from GNSS data. 
% Author :
% Jiangjun RAN, Pavel Ditmar
% Email:
% ranjj@sustech.edu.cn
% Note:
% Copyright (c) 2023 Jiangjun RAN. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% GPS_site_name_pavel={'KAGA';'MARG';'DKSG';'ASKY';'KULL';'SRMP';'RINK';'QAAR';'ILUL';'KELY';'KBUG';'NNVN';'TREO';'LYNS';'KULU';'KSNB';'PLPK';'SCOR';'GROK';'BLAS';'NRSK';'LEFN'};
% % apply without constraint
    GPS_site_name_pavel={'GMMA';'KAGA';'MARG';'DKSG';'ASKY';'KULL';'SRMP';'RINK';'QAAR';'ILUL';'KELY';'KBUG';'NNVN';'TREO';'LYNS';'KULU';'KSNB';'PLPK';'SCOR';'GROK';'BLAS';'NRSK';'LEFN';'QAQ1';'HEL2';'MIK2';'GMMA';'YMER';'LBIB';'WTHG';'HMBG';'DGJG';'VFDG';'HRDG';'KAPI';'SCBY';'KAGZ';'HEL2';'SENU';'HJOR';'KUAQ';'MIK2'};
    logic_constraint=0;
% % need constraint
%     GPS_site_name_pavel={'MARG'}; %'WTHG';'LBIB';'GROK';'KBUG';'BLAS';'DKSG';'KUAQ';'ASKY'
%     logic_constraint=1;
    
    number_scenario=0;
    
for kk=1:length(GPS_site_name_pavel)
    GPS_site_name=GPS_site_name_pavel{kk};
    load(['../Input/' GPS_site_name])
    E_r_group(1:number_scenario,1)=0.0;
    
    delta_tao=1.0/365.25
    n=length(day_decimal_year(:,1))
    for k=1:n-1
        liquid_p_time_series(k,2)=(liquid_p_time_series(k+1,2)-liquid_p_time_series(k,2))/delta_tao;
    end
    liquid_p_time_series(n,2)=0;
    for k=1:n-1
        RACMO_runoff_selected_period_interpolated(k,2)=(RACMO_runoff_selected_period_interpolated(k+1,2)-RACMO_runoff_selected_period_interpolated(k,2))/delta_tao;
    end
    RACMO_runoff_selected_period_interpolated(n,2)=0;
    
    figure
%     set(gcf,'WindowState','fullscreen')
    subplot(2,2,1)
    plot(day_decimal_year(:,1),RACMO_runoff_selected_period_interpolated(:,2),'k','linewidth',2)
    hold on
    %     axis([2009 2016 -100 10])
    grid on
    xlabel('Year')
    ylabel('Loading due to RACMO-based runoff [mm]')
    pbaspect([2 1 1])
    grid minor
    hold off
    
    if (isnan(GPS_time_series_running_mean_selected_period_interpolated_clean(1,2)))
        year_start=2010;
        year_end=2015; % up to 2016
        [GPS_time_series_running_mean_period_interpolated_clean_reduced]=selected_time_series_to_given_time_interval(year_start,year_end+1,GPS_time_series_running_mean_selected_period_interpolated_clean(:,1),GPS_time_series_running_mean_selected_period_interpolated_clean(:,2));
        day_decimal_year=[];
        day_decimal_year(:,1)=GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1);
        GPS_time_series_running_mean_selected_period_interpolated_clean=GPS_time_series_running_mean_period_interpolated_clean_reduced;
        [RACMO_SMB_selected_period_interpolated_reduced]=selected_time_series_to_given_time_interval(year_start,year_end+1,RACMO_SMB_selected_period_interpolated(:,1),RACMO_SMB_selected_period_interpolated(:,2));
        RACMO_SMB_selected_period_interpolated=RACMO_SMB_selected_period_interpolated_reduced;
        [liquid_p_time_series_reduced]=selected_time_series_to_given_time_interval(year_start,year_end+1,liquid_p_time_series(:,1),liquid_p_time_series(:,2));
        [RACMO_runoff_selected_period_interpolated_reduced]=selected_time_series_to_given_time_interval(year_start,year_end+1,RACMO_runoff_selected_period_interpolated(:,1),RACMO_runoff_selected_period_interpolated(:,2));
        
        if logic_constraint==1
%             [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_2022_3_11_ocean(2,2013,0.2,15,1.,GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1),GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
            [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig,formal_error_calibrated_E_r]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r(2,2013,0.,1.,15,1.,GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1),GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
            % peak mass per year
            meltwater_year=meltwater;
            meltwater_file='max_meltwater_year_';
%             meltwater_year=meltwater_orig; % RACMO meltwater production
%             meltwater_file='max_meltwater_orig_year_';
            max_meltwater_year=[];
            for k=year_start:year_end
                year_meltwater_re_selected=[];
                [year_meltwater_re_selected]=selected_time_series_to_given_time_interval(k,k+1,day_decimal_year(:,1),meltwater_year(:,1));
                [max_year,max_index] = max(abs(year_meltwater_re_selected(:,2)));
                max_meltwater_year(k-year_start+1,1)=year_meltwater_re_selected(max_index,1);
                max_meltwater_year(k-year_start+1,2)=year_meltwater_re_selected(max_index,2);
            end
            
            %% pavel uncertainty using Monte Carlo
            residual_GPS_SMB_fitted=GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2)-RACMO_SMB_selected_period_interpolated_reduced(:,2)-time_series_simulated;
            time_series_simulated_from_random_error=[];
%             invbeta_group(1:number_scenario,1)=0.0;
            invbeta_group(1:number_scenario,1)=0.0;
%             E_r_group(1:number_scenario,1)=0.0;
            for k=1:number_scenario
                random_error(:,1)=std(residual_GPS_SMB_fitted)*randn(1,length(residual_GPS_SMB_fitted));
                time_series_simulated_from_random_error(:,1)=GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1);
                time_series_simulated_from_random_error(:,2)=time_series_simulated+random_error+RACMO_SMB_selected_period_interpolated_reduced(:,2);
                [a_scenario,trend_scenario,c2_scenario,c3_scenario,E_r_scenario,theta_scenario,invbeta_scenario,time_series_simulated_scenario,discharge_scenario,meltwater_scenario,meltwater_orig_scenario,formal_error_calibrated_invbeta_scenario,calibrated_formal_error_meltwarer_scenario,intergated_upstream_runoff_per_year_scenario,intergated_upstream_runoff_per_year_orig_scenario]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean(2,2013,0.,1.,15,1.,time_series_simulated_from_random_error(:,1),time_series_simulated_from_random_error(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
                invbeta_group(k,1)=invbeta_scenario;
                
                %% scale factor Monte Carlo
                if year_start==2009
                    E_r_group(k,1)=1.0+E_r_scenario(2,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(3,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(4,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(5,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(6,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(7,1); %2015
                elseif year_start==2010
                    E_r_group(k,1)=1.0+E_r_scenario(1,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(2,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(3,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(4,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(5,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(6,1); %2015
                end
            end
        else
%             [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_2022_3_11_ocean(2,0.,0.2,15,1.,GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1),GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
            [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig,formal_error_calibrated_E_r]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r(2,0.,0.,1.,15,1.,GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1),GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
            % peak mass per year
            meltwater_year=meltwater;
            meltwater_file='max_meltwater_year_';
%             meltwater_year=meltwater_orig; % RACMO meltwater production
%             meltwater_file='max_meltwater_orig_year_';
            max_meltwater_year=[];
            for k=year_start:year_end
                year_meltwater_re_selected=[];
                [year_meltwater_re_selected]=selected_time_series_to_given_time_interval(k,k+1,day_decimal_year(:,1),meltwater_year(:,1));
                [max_year,max_index] = max(abs(year_meltwater_re_selected(:,2)));
                max_meltwater_year(k-year_start+1,1)=year_meltwater_re_selected(max_index,1);
                max_meltwater_year(k-year_start+1,2)=year_meltwater_re_selected(max_index,2);
            end
            
            %% pavel uncertainty using Monte Carlo
            residual_GPS_SMB_fitted=GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2)-RACMO_SMB_selected_period_interpolated_reduced(:,2)-time_series_simulated;
            time_series_simulated_from_random_error=[];
%             invbeta_group(1:number_scenario,1)=0.0;
            invbeta_group(1:number_scenario,1)=0.0;
            for k=1:number_scenario
                random_error(:,1)=std(residual_GPS_SMB_fitted)*randn(1,length(residual_GPS_SMB_fitted));
                time_series_simulated_from_random_error(:,1)=GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1);
                time_series_simulated_from_random_error(:,2)=time_series_simulated+random_error+RACMO_SMB_selected_period_interpolated_reduced(:,2);
                [a_scenario,trend_scenario,c2_scenario,c3_scenario,E_r_scenario,theta_scenario,invbeta_scenario,time_series_simulated_scenario,discharge_scenario,meltwater_scenario,meltwater_orig_scenario,formal_error_calibrated_invbeta_scenario,calibrated_formal_error_meltwarer_scenario,intergated_upstream_runoff_per_year_scenario,intergated_upstream_runoff_per_year_orig_scenario]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean(2,0,0.,1.,15,1.,time_series_simulated_from_random_error(:,1),time_series_simulated_from_random_error(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
                invbeta_group(k,1)=invbeta_scenario;
                
                %% scale factor Monte Carlo
                if year_start==2009
                    E_r_group(k,1)=1.0+E_r_scenario(2,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(3,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(4,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(5,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(6,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(7,1); %2015
                elseif year_start==2010
                    E_r_group(k,1)=1.0+E_r_scenario(1,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(2,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(3,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(4,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(5,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(6,1); %2015
                end
            end
        end
    else
        year_start=2009;
        year_end=2015; % up to 2016
        if logic_constraint==1
%            [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_2022_3_11_ocean(2,2013,0.,1.,15,1.,GPS_time_series_running_mean_period_interpolated_clean_reduced(:,1),GPS_time_series_running_mean_period_interpolated_clean_reduced(:,2),RACMO_SMB_selected_period_interpolated_reduced(:,2),liquid_p_time_series_reduced(:,2),RACMO_runoff_selected_period_interpolated_reduced(:,2),year_start,year_end);
            [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig,formal_error_calibrated_E_r]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r(2,2013,0.,1.,15,1.,GPS_time_series_running_mean_selected_period_interpolated_clean(:,1),GPS_time_series_running_mean_selected_period_interpolated_clean(:,2),RACMO_SMB_selected_period_interpolated(:,2),liquid_p_time_series(:,2),RACMO_runoff_selected_period_interpolated(:,2),year_start,year_end);
            % peak mass per year
            meltwater_year=meltwater;
            meltwater_file='max_meltwater_year_';
%             meltwater_year=meltwater_orig; % RACMO meltwater production
%             meltwater_file='max_meltwater_orig_year_';
            max_meltwater_year=[];
            for k=year_start:year_end
                year_meltwater_re_selected=[];
                [year_meltwater_re_selected]=selected_time_series_to_given_time_interval(k,k+1,day_decimal_year(:,1),meltwater_year(:,1));
                [max_year,max_index] = max(abs(year_meltwater_re_selected(:,2)));
                max_meltwater_year(k-year_start+1,1)=year_meltwater_re_selected(max_index,1);
                max_meltwater_year(k-year_start+1,2)=year_meltwater_re_selected(max_index,2);
            end
            
            %% pavel uncertainty using Monte Carlo
            residual_GPS_SMB_fitted=GPS_time_series_running_mean_selected_period_interpolated_clean(:,2)-RACMO_SMB_selected_period_interpolated(:,2)-time_series_simulated;
            time_series_simulated_from_random_error=[];
%             invbeta_group(1:number_scenario,1)=0.0;
            invbeta_group(1:number_scenario,1)=0.0;
            for k=1:number_scenario
                random_error(:,1)=std(residual_GPS_SMB_fitted)*randn(1,length(residual_GPS_SMB_fitted));
                time_series_simulated_from_random_error(:,1)=GPS_time_series_running_mean_selected_period_interpolated_clean(:,1);
                time_series_simulated_from_random_error(:,2)=time_series_simulated+random_error+RACMO_SMB_selected_period_interpolated(:,2);
                [a_scenario,trend_scenario,c2_scenario,c3_scenario,E_r_scenario,theta_scenario,invbeta_scenario,time_series_simulated_scenario,discharge_scenario,meltwater_scenario,meltwater_orig_scenario,formal_error_calibrated_invbeta_scenario,calibrated_formal_error_meltwarer_scenario,intergated_upstream_runoff_per_year_scenario,intergated_upstream_runoff_per_year_orig_scenario]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean(2,2013,0.,1.,15,1.,time_series_simulated_from_random_error(:,1),time_series_simulated_from_random_error(:,2),RACMO_SMB_selected_period_interpolated(:,2),liquid_p_time_series(:,2),RACMO_runoff_selected_period_interpolated(:,2),year_start,year_end);
                invbeta_group(k,1)=invbeta_scenario;
                
                %% scale factor Monte Carlo
                if year_start==2009
                    E_r_group(k,1)=1.0+E_r_scenario(2,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(3,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(4,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(5,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(6,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(7,1); %2015
                elseif year_start==2010
                    E_r_group(k,1)=1.0+E_r_scenario(1,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(2,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(3,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(4,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(5,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(6,1); %2015
                end
            end
        else
%             [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_2022_3_11_ocean(2,0.,0.2,15,1.,day_decimal_year(:,1),GPS_time_series_running_mean_selected_period_interpolated_clean(:,2),RACMO_SMB_selected_period_interpolated(:,2),liquid_p_time_series(:,2),RACMO_runoff_selected_period_interpolated(:,2),year_start,year_end);
            [a,trend,c2,c3,E_r,theta,invbeta,time_series_simulated,discharge,meltwater,meltwater_orig,formal_error_calibrated_invbeta,calibrated_formal_error_meltwarer,intergated_upstream_runoff_per_year,intergated_upstream_runoff_per_year_orig,formal_error_calibrated_E_r]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean_E_r(2,0.,0.,1.,15,1.,GPS_time_series_running_mean_selected_period_interpolated_clean(:,1),GPS_time_series_running_mean_selected_period_interpolated_clean(:,2),RACMO_SMB_selected_period_interpolated(:,2),liquid_p_time_series(:,2),RACMO_runoff_selected_period_interpolated(:,2),year_start,year_end);
            % peak mass per year
            meltwater_year=meltwater;
            meltwater_file='max_meltwater_year_';
%             meltwater_year=meltwater_orig; % RACMO meltwater production
%             meltwater_file='max_meltwater_orig_year_';
            max_meltwater_year=[];
            for k=year_start:year_end
                year_meltwater_re_selected=[];
                [year_meltwater_re_selected]=selected_time_series_to_given_time_interval(k,k+1,day_decimal_year(:,1),meltwater_year(:,1));
                [max_year,max_index] = max(abs(year_meltwater_re_selected(:,2)));
                max_meltwater_year(k-year_start+1,1)=year_meltwater_re_selected(max_index,1);
                max_meltwater_year(k-year_start+1,2)=year_meltwater_re_selected(max_index,2);
            end
            
            %% pavel uncertainty using Monte Carlo
            residual_GPS_SMB_fitted=GPS_time_series_running_mean_selected_period_interpolated_clean(:,2)-RACMO_SMB_selected_period_interpolated(:,2)-time_series_simulated;
            time_series_simulated_from_random_error=[];
%             invbeta_group(1:number_scenario,1)=0.0;
            invbeta_group(1:number_scenario,1)=0.0;
            for k=1:number_scenario
                random_error(:,1)=std(residual_GPS_SMB_fitted)*randn(1,length(residual_GPS_SMB_fitted));
                time_series_simulated_from_random_error(:,1)=GPS_time_series_running_mean_selected_period_interpolated_clean(:,1);
                time_series_simulated_from_random_error(:,2)=time_series_simulated+random_error+RACMO_SMB_selected_period_interpolated(:,2);
                [a_scenario,trend_scenario,c2_scenario,c3_scenario,E_r_scenario,theta_scenario,invbeta_scenario,time_series_simulated_scenario,discharge_scenario,meltwater_scenario,meltwater_orig_scenario,formal_error_calibrated_invbeta_scenario,calibrated_formal_error_meltwarer_scenario,intergated_upstream_runoff_per_year_scenario,intergated_upstream_runoff_per_year_orig_scenario]=est_Thet_Bet_y1Er_Reg0Eanyy_Reg0B_Reg0ErMean(2,0,0.,1.,15,1.,time_series_simulated_from_random_error(:,1),time_series_simulated_from_random_error(:,2),RACMO_SMB_selected_period_interpolated(:,2),liquid_p_time_series(:,2),RACMO_runoff_selected_period_interpolated(:,2),year_start,year_end);
                invbeta_group(k,1)=invbeta_scenario;
                
                %% scale factor Monte Carlo
                if year_start==2009
                    E_r_group(k,1)=1.0+E_r_scenario(2,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(3,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(4,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(5,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(6,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(7,1); %2015
                elseif year_start==2010
                    E_r_group(k,1)=1.0+E_r_scenario(1,1); %2010
                    E_r_group(k,2)=1.0+E_r_scenario(2,1); %2011
                    E_r_group(k,3)=1.0+E_r_scenario(3,1); %2012
                    E_r_group(k,4)=1.0+E_r_scenario(4,1); %2013
                    E_r_group(k,5)=1.0+E_r_scenario(5,1); %2014
                    E_r_group(k,6)=1.0+E_r_scenario(6,1); %2015
                end
            end
        end
    end
    
    %% scale factor to RACMO runoff
    if year_start==2009
        scale_factor_to_runoff{kk,1}=GPS_site_name;
        scale_factor_to_runoff{kk,2}=1.0+E_r(2,1); %2010
        scale_factor_to_runoff{kk,3}=1.0+E_r(3,1); %2011
        scale_factor_to_runoff{kk,4}=1.0+E_r(4,1); %2012
        scale_factor_to_runoff{kk,5}=1.0+E_r(5,1); %2013
        scale_factor_to_runoff{kk,6}=1.0+E_r(6,1); %2014
        scale_factor_to_runoff{kk,7}=1.0+E_r(7,1); %2015
        
        scale_factor_to_runoff_error{kk,1}=GPS_site_name;
        scale_factor_to_runoff_error{kk,2}=formal_error_calibrated_E_r(2,1); %2010
        scale_factor_to_runoff_error{kk,3}=formal_error_calibrated_E_r(3,1); %2011
        scale_factor_to_runoff_error{kk,4}=formal_error_calibrated_E_r(4,1); %2012
        scale_factor_to_runoff_error{kk,5}=formal_error_calibrated_E_r(5,1); %2013
        scale_factor_to_runoff_error{kk,6}=formal_error_calibrated_E_r(6,1); %2014
        scale_factor_to_runoff_error{kk,7}=formal_error_calibrated_E_r(7,1); %2015
    elseif year_start==2010
        scale_factor_to_runoff{kk,1}=GPS_site_name;
        scale_factor_to_runoff{kk,2}=1.0+E_r(1,1); %2010
        scale_factor_to_runoff{kk,3}=1.0+E_r(2,1); %2011
        scale_factor_to_runoff{kk,4}=1.0+E_r(3,1); %2012
        scale_factor_to_runoff{kk,5}=1.0+E_r(4,1); %2013
        scale_factor_to_runoff{kk,6}=1.0+E_r(5,1); %2014
        scale_factor_to_runoff{kk,7}=1.0+E_r(6,1); %2015
        
        scale_factor_to_runoff_error{kk,1}=GPS_site_name;
        scale_factor_to_runoff_error{kk,2}=formal_error_calibrated_E_r(1,1); %2010
        scale_factor_to_runoff_error{kk,3}=formal_error_calibrated_E_r(2,1); %2011
        scale_factor_to_runoff_error{kk,4}=formal_error_calibrated_E_r(3,1); %2012
        scale_factor_to_runoff_error{kk,5}=formal_error_calibrated_E_r(4,1); %2013
        scale_factor_to_runoff_error{kk,6}=formal_error_calibrated_E_r(5,1); %2014
        scale_factor_to_runoff_error{kk,7}=formal_error_calibrated_E_r(6,1); %2015
    end
    
%     figure
    subplot(2,2,2)
    plot(day_decimal_year(:,1),GPS_time_series_running_mean_selected_period_interpolated_clean(:,2)-RACMO_SMB_selected_period_interpolated(:,2),day_decimal_year(:,1),time_series_simulated,'k','linewidth',2)
    hold on
    %     axis([2009 2016 -50 50])
    grid on
    xlabel('Year')
    ylabel('Residual loading: observed vs simulated [mm]')
    pbaspect([2 1 1])
    grid minor
    hold off
    
%     figure
    subplot(2,2,3)
    plot(day_decimal_year(:,1),discharge,'k','linewidth',2)
    title(['invbeta' num2str(invbeta*365.25) 'error' num2str(formal_error_calibrated_invbeta*365.25)])
    hold on

    grid on
    xlabel('Year')
    ylabel('Discharge-related loading [mm]')
    pbaspect([2 1 1])
    grid minor
    hold off
    
%     figure
    subplot(2,2,4)
    plot(day_decimal_year(:,1),meltwater_orig,'r',day_decimal_year(:,1),meltwater,'k','linewidth',2)
    hold on
    plot(max_meltwater_year(:,1),max_meltwater_year(:,2),'g-p')
    xlim([2009 2016])
    grid on
    xlabel('Year')
    ylabel('Liquid water loading [mm]')
    pbaspect([2 1 1])
    grid minor
    hold off
    
    title(GPS_site_name)
    if logic_constraint==1
        file_in = ['..\Output\' 'constraint\meltwater_loading_modified_' GPS_site_name '.png'];
        saveas(gcf,file_in,'png')
        file_in = ['..\Output\' 'constraint\meltwater_loading_modefied_' GPS_site_name '.pdf'];
        saveas(gcf,file_in,'pdf')
    else
        file_in = ['..\Output\' 'meltwater_loading_modified_' GPS_site_name '.png'];
        saveas(gcf,file_in,'png')
        file_in = ['..\Output\' 'meltwater_loading_modefied_' GPS_site_name '.pdf'];
        saveas(gcf,file_in,'pdf')
    end
    
    figure
	year_plotted=[1:year_end-year_start+1];
	plot([year_start:1:year_end],intergated_upstream_runoff_per_year_orig(end,year_plotted),'r',[year_start:1:year_end],intergated_upstream_runoff_per_year(end,year_plotted),'k','linewidth',2)
	legend('intergated\_upstream\_runoff\_orig','intergated\_upstream\_runoff')
	grid on
	xlabel('Year')
	ylabel_text=['Loading due to integrated upstream runoff in year ',int2str(year_plotted),' [mm]'];
	ylabel(ylabel_text)
	pbaspect([2 1 1])
	grid minor
	hold off
    title(GPS_site_name)
    file_in = ['..\Output\' 'Loading due to integrated upstream runoff in year ' int2str(year_plotted) GPS_site_name '.png'];
    saveas(gcf,file_in,'png')
    file_in = ['..\Output\' 'Loading due to integrated upstream runoff in year ' int2str(year_plotted) GPS_site_name '.pdf'];
    saveas(gcf,file_in,'pdf')
    file_simulated_mass = ['..\Output\' 'Loading due to integrated upstream runoff_' GPS_site_name '.txt'];
    fid=fopen(file_simulated_mass,'w');
    for k=year_start:year_end
        fprintf(fid,'%u     %u      %u      %u\r\n',k,intergated_upstream_runoff_per_year(end,k-year_start+1),intergated_upstream_runoff_per_year_orig(end,k-year_start+1),intergated_upstream_runoff_per_year(end,k-year_start+1)-intergated_upstream_runoff_per_year_orig(end,k-year_start+1));
    end
    fclose(fid);
    
    figure
    calibrated_formal_error_meltwarer(1:1:length(day_decimal_year),1)=mean(calibrated_formal_error_meltwarer(1:1:length(day_decimal_year),1));
    %errorbar(day_decimal_year(1:1:length(day_decimal_year),1),meltwater(1:1:length(day_decimal_year),1),calibrated_formal_error_meltwarer(1:1:length(day_decimal_year),1),'r')
    h8=patch([day_decimal_year(:,1)' day_decimal_year(length(day_decimal_year):-1:1,1)'],[meltwater(:,1)'+3*calibrated_formal_error_meltwarer(:,1)' meltwater(length(day_decimal_year):-1:1,1)'-3*calibrated_formal_error_meltwarer(length(day_decimal_year):-1:1,1)'],[0.7 0.7 0.7],'EdgeColor', [0.7 0.7 0.7]);
    set(h8,'FaceColor',[1 0.800000011920929 0.800000011920929],'FaceAlpha',0.8)
    hold on
    plot(day_decimal_year(:,1),meltwater(:,1),'r','Linewidth',1.5)
    xlim([2009 2016])
    grid on
    xlabel('Year')
    ylabel('Liquid water loading [mm]')
    pbaspect([3 1 1])
    grid minor
    hold off
    title(GPS_site_name)
    
    % Scale file
    file_error = ['..\Output\' 'Scale_file' num2str(number_scenario) '2010_2015' GPS_site_name 'error'  '.txt'];
    fid=fopen(file_error,'w');
    fprintf(fid,'%s      \r\n','2010 2011 2012 2013 2014 2015');
    fprintf(fid,'%u     %u      %u      %u      %u      %u      \r\n',scale_factor_to_runoff{kk,2},scale_factor_to_runoff{kk,3},scale_factor_to_runoff{kk,4},scale_factor_to_runoff{kk,5},scale_factor_to_runoff{kk,6},scale_factor_to_runoff{kk,7});
    fprintf(fid,'%u     %u      %u      %u      %u      %u      \r\n',scale_factor_to_runoff_error{kk,2},scale_factor_to_runoff_error{kk,3},scale_factor_to_runoff_error{kk,4},scale_factor_to_runoff_error{kk,5},scale_factor_to_runoff_error{kk,6},scale_factor_to_runoff_error{kk,7});
    fclose(fid);
    % Monte Carlo scenarios
    if number_scenario>0
        figure
        hist(invbeta_group*365.25-invbeta*365.25)
        xlabel('invbeta_group-invbeta (day)')
        MonteCarlo_residual=sort(invbeta_group*365.25-invbeta*365.25,'ascend');
        title(['invbeta MonteCarlo_residual' num2str(invbeta*365.25) 'error' num2str(MonteCarlo_residual(26)) 'GRL' num2str(MonteCarlo_residual(76))])
        file_in = ['..\Output\' 'meltwater storage uncertainty '  GPS_site_name '.png'];
        saveas(gcf,file_in,'png')
        file_in = ['..\Output\' 'meltwater storage uncertainty ' GPS_site_name '.pdf'];
        saveas(gcf,file_in,'pdf')
        
        % Scale file
        file_scenario = ['..\Output\' 'Scale_file' num2str(number_scenario) 'scanarios' '2010_2015' GPS_site_name '.txt'];
        fid=fopen(file_scenario,'w');
        fprintf(fid,'%s      \r\n','2010 2011 2012 2013 2014 2015');
        fprintf(fid,'%u     %u      %u      %u      %u      %u      \r\n',scale_factor_to_runoff{kk,2},scale_factor_to_runoff{kk,3},scale_factor_to_runoff{kk,4},scale_factor_to_runoff{kk,5},scale_factor_to_runoff{kk,6},scale_factor_to_runoff{kk,7});
        for k=1:number_scenario
            fprintf(fid,'%u     %u      %u      %u      %u      %u      \r\n',E_r_group(k,1),E_r_group(k,2),E_r_group(k,3),E_r_group(k,4),E_r_group(k,5),E_r_group(k,6));
        end
        fclose(fid);
    end
    
    if logic_constraint==1
        file_in = ['..\Output\' 'constraint\meltwater_loading_modified_error_bar_zone_' GPS_site_name '.png'];
        saveas(gcf,file_in,'png')
        file_in = ['..\Output\' 'constraint\meltwater_loading_modefied_error_bar_zone_' GPS_site_name '.pdf'];
        saveas(gcf,file_in,'pdf')
        file_simulated_mass = ['meltwater_loading_pavel\' 'constraint\meltwater_loading_modified_error_bar_zone_' GPS_site_name '.txt'];
        fid=fopen(file_simulated_mass,'w');
        for k=1:length(meltwater(:,1))
            fprintf(fid,'%u     %u      %u\r\n',day_decimal_year(k,1),meltwater(k,1),calibrated_formal_error_meltwarer(k,1));
        end
        fclose(fid);
        file_simulated_mass = ['..\Output\' meltwater_file GPS_site_name '.txt'];
        fid=fopen(file_simulated_mass,'w');
        for k=1:length(max_meltwater_year(:,1))
            fprintf(fid,'%u     %u      \r\n',max_meltwater_year(k,1),max_meltwater_year(k,2));
        end
        fclose(fid);
    else
        file_in = ['..\Output\' 'meltwater_loading_modified_error_bar_zone_' GPS_site_name '.png'];
        saveas(gcf,file_in,'png')
        file_in = ['..\Output\' 'meltwater_loading_modefied_error_bar_zone_' GPS_site_name '.pdf'];
        saveas(gcf,file_in,'pdf')
        file_simulated_mass = ['..\Output\' 'meltwater_loading_modified_error_bar_zone_' GPS_site_name '.txt'];
        fid=fopen(file_simulated_mass,'w');
        for k=1:length(meltwater(:,1))
            fprintf(fid,'%u     %u      %u\r\n',day_decimal_year(k,1),meltwater(k,1),calibrated_formal_error_meltwarer(k,1));
        end
        fclose(fid);
        file_simulated_mass = ['..\Output\' meltwater_file GPS_site_name '.txt'];
        fid=fopen(file_simulated_mass,'w');
        for k=1:length(max_meltwater_year(:,1))
            fprintf(fid,'%u     %u      \r\n',max_meltwater_year(k,1),max_meltwater_year(k,2));
        end
        fclose(fid);
    end
    
    close all
    clearvars -except GPS_site_name_pavel kk logic_constraint scale_factor_to_runoff scale_factor_to_runoff_error number_scenario
end