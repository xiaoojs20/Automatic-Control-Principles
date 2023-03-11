function sys=nyquist_freqscan(sys,state,w)
    freqresp_sys = reshape(freqresp(sys,w),[size(w),1]);
    abs_sys = abs(freqresp_sys);
    angle_sys = angle(freqresp_sys);
    z_sys=abs_sys.*exp(1i.*angle_sys);
    
    figure;
    plot(z_sys,'b.-','LineWidth',2);
    grid on;
    xlabel('Re','FontSize',12,'Interpreter','Latex');
    ylabel('Im','FontSize',12,'Interpreter','Latex');
    
    if state == "delta"
        title("The Nyquist diagram of $G_{\delta}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    elseif state == "omega"
        title("The Nyquist diagram of $G_{\omega}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    elseif state == "Pem"
        title("The Nyquist diagram of $G_{P_{em}}(j\omega)$ (frequency scan)",'FontSize',16,'Interpreter','latex');
    end
    sys=sys;

end