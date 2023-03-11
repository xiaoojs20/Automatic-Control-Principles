function sys=bode_freqscan(sys,state,w)
    freqresp_sys = reshape(freqresp(sys,w),[size(w),1]);
    abs_sys = abs(freqresp_sys);
    angle_sys = angle(freqresp_sys);
    
    figure;
    subplot(2,1,1);
    semilogx(w,20.*log10(abs_sys),'b.-','LineWidth',2);
    grid on;
    xlabel('$\omega$ /(rad/s)','FontSize',12,'Interpreter','Latex');
    ylabel('$L(\omega)$/dB','FontSize',12,'Interpreter','Latex');
    xlim([10^-2,10^2]);
    subplot(2,1,2);
    semilogx(w,angle_sys.*180./pi,'b.-','LineWidth',2);
    xlim([10^-2,10^2]);
    grid on;
    xlabel('$\omega$ /(rad/s)','FontSize',12,'Interpreter','Latex');
    ylabel('$\phi(\omega)/^{\circ}$','FontSize',12,'Interpreter','Latex');

    if state == "delta"
        sgtitle("The Bode diagram of $G_{\delta}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    elseif state == "omega"
        sgtitle("The Bode diagram of $G_{\omega}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    elseif state == "Pem"
        sgtitle("The Bode diagram of $G_{P_{em}}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    end
    sys=sys;

end