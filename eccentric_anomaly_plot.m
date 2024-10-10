function eccentric_anomaly_plot(a,e,mu,N,k,E0,t0)
    T = 2*pi*sqrt(a^3/mu);
    t = linspace(t0,k*T,N);
    legend_names = strings(0);

    for i = 1 : length(e)
        E = kepler_equation(t,e(i),a,mu,t0,E0);
        plot(t./T,rad2deg(E),'LineWidth',2); hold on; grid on;
        legend_names(i) = ['e = ' num2str(e(i))];
    end
    legend(legend_names);
    xlim([0 t(end)./T]); ylim([0 max(rad2deg(E))]); yticks([linspace(0,max(rad2deg(E)),7)])
    title('Eccentric anomaly - Time'); xlabel('t [T]'); ylabel('E [deg]');
end