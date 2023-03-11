function sys=abs_angle_freqscan(sys,state,w)
    freqresp_sys = reshape(freqresp(sys,w),[size(w),1]);
    abs_sys = abs(freqresp_sys);
    angle_sys = angle(freqresp_sys);
    if size(w)==2000
        w_axis=[-1000:-1, 1:1000];
    else
        w_axis= 1:1000;
    end
%     figure;
    plot(w_axis,abs_sys,'LineWidth',2);
    hold on;
    grid on;
    xlabel('$\omega$ /(rad/s)','FontSize',12,'Interpreter','Latex');
    if state == "delta"
        ylabel('$\left| G_{\delta}(j\omega) \right|$','FontSize',12,'Interpreter','Latex');
        title("The amplitude of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
    elseif state == "omega"
        ylabel('$\left| G_{\omega}(j\omega) \right|$','FontSize',12,'Interpreter','Latex');
        title("The amplitude of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
    elseif state == "Pem"
        ylabel('$\left| G_{P_{em}}(j\omega) \right|$','FontSize',12,'Interpreter','Latex');
        title("The amplitude of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');
    end

%     figure;
%     plot(w_axis,angle_sys.*180./pi,'b.-','LineWidth',2);
%     grid on;
%     xlabel('$\omega$ /(rad/s)','FontSize',12,'Interpreter','Latex');
%     if state == "delta"
%         ylabel('$\angle G_{\delta}(j\omega)$','FontSize',12,'Interpreter','Latex');
%         title("The angle of $G_{\delta}(j\omega)$",'FontSize',16,'Interpreter','latex');
%     elseif state == "omega"
%         ylabel('$\angle G_{\omega}(j\omega)$','FontSize',12,'Interpreter','Latex');
%         title("The angle of $G_{\omega}(j\omega)$",'FontSize',16,'Interpreter','latex');
%     elseif state == "Pem"
%         ylabel('$\angle G_{P_{em}}(j\omega)$','FontSize',12,'Interpreter','Latex');
%         title("The angle of $G_{P_{em}}(j\omega)$",'FontSize',16,'Interpreter','latex');
%     end

    sys=sys;

end