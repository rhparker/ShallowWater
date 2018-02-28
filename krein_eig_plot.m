% plot Krein eigenvalues for values in R+

% min_x = 0.0047731;
% max_x = 0.0047733;
% num_pts = 50;

min_x = 15.80;
max_x = 16.00;
num_pts = 100;

x = linspace(min_x, max_x, num_pts);

% % figure out where determinant of Krein matrix is 0
% z = fsolve( @(z) det(kreinmatrix_z(z)), 0.0691^2) 

kEigs = [];
for index = 1:num_pts
    [Kz, lambda] = kreinmatrix_z( x(index) );
    kEigs = [kEigs lambda];
    disp(index);
end

if length(lambda) == 2
    bestfit1 = fit(x', kEigs(1,:)', 'poly1');
    bestfit2 = fit(x', kEigs(2,:)', 'poly1');

    figure;
    hold on;
    plot(x, kEigs(1,:), '.r', x, kEigs(2,:), '.b')
%     plot(x, bestfit1(x), 'r');
%     plot(x, bestfit2(x), 'b');
    title('Krein eigenvalues, Double Pulse 2(3)');
    legend('Krein Eigenvalue 1', 'Krein Eigenvalue 2');
    axis([15.8 16 -5 5]);

    figure;
    hold on;
    plot(x, kEigs(2,:), '.b')
%     plot(x, bestfit2(x), 'b');
    title('Krein eigenvalue 2, Double Pulse 2(3)');
    legend('Krein Eigenvalue 2');
    axis([15.8 16 -2e-3 2e-3]);
    
%     % plot with different y axes
%     figure;
%     [hAx,hLine1,hLine2] = plotyy(x, kEigs(1,:), x, kEigs(2,:) );
%     set(hAx(1),'ycolor','r') 
%     set(hAx(2),'ycolor','b')
%     set(hLine1,'marker','.','color','r','linestyle','none')
%     set(hLine2,'marker','.','color','b','linestyle','none')
%     maxval = cellfun(@(x) max(abs(x)), get([hLine1 hLine2], 'YData'));
%     ylim = [-maxval, maxval] * 1.1;  % Mult by 1.1 to pad out a bit
%     set(hAx(1), 'YLim', ylim(1,:) / 50 );
%     set(hAx(2), 'YLim', ylim(2,:) );
%     title('Krein eigenvalue, different y axis scales, Double Pulse 2(3)');
%     legend('Krein Eigenvalue 1', 'Krein Eigenvalue 2');
    
    
    P = InterX([x;kEigs(1,:)],[x;kEigs(2,:)]); 
    PZ1 = InterX([x;kEigs(1,:)],[x;0*x]);
    PZ2 = InterX([x;kEigs(2,:)],[x;0*x]);
    
    
else
    bestfit1 = fit(x', kEigs(1,:)', 'poly1');

    figure;
    hold on;
    plot(x, kEigs(1,:), '.r');
    plot(x, bestfit1(x), 'r');
    title('Krein eigenvalue, Double Pulse 2(2)');
    legend('Krein Eigenvalue');
    
    P = InterX([x;kEigs(1,:)],[x;0*x]); 
end
