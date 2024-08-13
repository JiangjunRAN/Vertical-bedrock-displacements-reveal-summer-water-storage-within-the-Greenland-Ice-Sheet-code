function [time_series_re_selected]=selected_time_series_to_given_time_interval(time_strat,time_end,time,data)
    n=1;
    for k=1:length(time)
        if (time(k)>=time_strat && time(k)<=time_end)
            time_series_re_selected(n,1)=time(k);
            time_series_re_selected(n,2)=data(k);
            n=n+1;
        end
    end
end