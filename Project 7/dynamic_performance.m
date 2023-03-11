%% 
% parameters = [overshoot damping_factor stable_time]
% This function, with input time and system response values, 
% outputs the dynamic characteristics of the system, such as
% overshoot, damping_factor, stable_time.
% 
% Xiao Jinsong
% 2022/10/25

%% Function
function [overshoot, damping_factor, stable_time]=dynamic_performance(ts,values,mode)
    if mode==0
        model_name = "linear model";
    end
    if mode==1
        model_name = "nonlinear model";
    end
    tolerance = 0.02;
    figure;
    plot(ts,values,'LineWidth',2);
    steady_period = 20;
    steady_period_max = max(values(length(ts)-steady_period:length(ts)));
    steady_period_min = min(values(length(ts)-steady_period:length(ts)));
    steadystate_value = (steady_period_max+steady_period_min)/2;
    peak_value = max(values);
    uppper_value = steadystate_value*(1+tolerance);
    lower_value = steadystate_value*(1-tolerance);
    peak_values = repmat(peak_value,length(ts),1);
    steadystate_values = repmat(steadystate_value,length(ts),1);
    uppper_values = repmat(uppper_value,length(ts),1);
    lower_values = repmat(lower_value,length(ts),1);
    hold on;
    grid on;
    plot(ts,steadystate_values,'LineWidth',1,'Color','r');
    plot(ts,peak_values,'LineWidth',1,'LineStyle','--','Color','g')
    plot(ts,uppper_values,'LineWidth',1,'LineStyle','--','Color','c')
    plot(ts,lower_values,'LineWidth',1,'LineStyle','--','Color','m')
    hold off;
    
    xlabel('time/s','FontSize',16,'Interpreter','Latex');
    ylabel('$V_{t}/V$','FontSize',16,'Interpreter','Latex');
    title(['Dynamic performance analysis of ' model_name],'FontSize',16,'Interpreter','Latex');
    legend("$V_{t}$", "steady state value", "peak value", "upper value", "lower value",...
        'Interpreter','latex','Location','SouthEast');
    overshoot = (peak_value - steadystate_value)/steadystate_value;
    
    for stime = 1:0.1:100   
        after_stime = ts(ts>stime);
        after_stime_num = find(ts==after_stime(1));
        after_stime_value = values(after_stime_num:end);
        after_time_max = max(after_stime_value);
        after_time_min = min(after_stime_value);
        if((after_time_max<=uppper_value) && after_time_min>=lower_value)
            stable_time = stime;
            break;
        end
    end
    damping_factor = sqrt((log(overshoot)^2)/(pi^2+log(overshoot)^2));
 
    fprintf("~~~~~~~> The overshoot of %s is %f percent. <~~~~~~~\n",model_name,overshoot*100);
    fprintf("~~~~~~~> The damping factor of %s is %f. <~~~~~~~\n",model_name,damping_factor);
    fprintf("~~~~~~~> The stable time of %s is %f. <~~~~~~~\n",model_name,stable_time);

end