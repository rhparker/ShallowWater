% use best-fit line to approximate speed
% of double pulse complex

function speed = center_speed(time, centers, cutoff)

    if ~exist('cutoff','var')
        cutoff = length(time);
    end
    
    figure;
    hold on;
    xplot = time(1:cutoff);
    yplot = centers(1:cutoff);
    
    marker_size = 15;
    scatter(xplot, yplot, marker_size, 'filled');

    % best fit line 
    bestfit = fit(xplot', yplot', 'poly1');
    plot(xplot, bestfit(xplot));
    plot_title = 'location of center of double pulse complex vs time';
    title({plot_title, strcat('slope =  ',num2str(abs(bestfit.p1)), '   intercept =  ',num2str(abs(bestfit.p2)))});
    ylabel('center of double pulse complex');
    xlabel('time');

    speed = bestfit.p1
end