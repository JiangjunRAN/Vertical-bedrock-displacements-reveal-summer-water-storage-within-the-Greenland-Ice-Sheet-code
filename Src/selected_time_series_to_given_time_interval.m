function [time_series_re_selected]=selected_time_series_to_given_time_interval(time_strat,time_end,time,data)
%***********************************************************************************************%
% Copyright (c) 2024-2099 by Gravity and Cryosphere Lab @ SUSTecch, all rights reserved.
% Author: Jiangjun Ran
% E-mail: ranjj@sustech.edu.cn
% cite: Ran, J., Ditmar, P., van den Broeke, M. R., et al. (in press). Vertical bedrock shifts reveal summer water storage in Greenland ice sheet. Nature.
%***********************************************************************************************%
    n=1;
    for k=1:length(time)
        if (time(k)>=time_strat && time(k)<=time_end)
            time_series_re_selected(n,1)=time(k);
            time_series_re_selected(n,2)=data(k);
            n=n+1;
        end
    end
end