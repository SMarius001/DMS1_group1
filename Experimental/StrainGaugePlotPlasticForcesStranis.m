% clear; clc; close all;
% 
% T = readtable("../DataComparisonStrainGaugesandForce.xlsx",Range="FA3:FI9074");
% 
% Adr1 = ["R0","R90","R45";"A0","A90","";"S2","S1",""];
% Adr2 = [3, 2, 2];
% 
% for i = 1:3
%     for j = 1:Adr2(i)
%         Temp = T{:,Adr1(i,j)};
%         k = 1;
%         while abs(Temp(k)) < 19300
%             k = k + 1;
%         end
%         T{k+1:end,Adr1(i,j)} = NaN;
%     end
% end
% 
% for i = 1:3
%     figure('Color','w')
%     hold on
%     for j = 1:Adr2(i)
%         plot(T,Adr1(i,j),"F")
%     end
%     hold off
%     xlabel("Microstrain")
%     ylabel("Force [kN]")
%     xtickformat('%,g')
%     xticks(-20000:10000:20000)
%     ax = gca; 
%     currentTicks = get(ax, 'XTick');
%     labels = strtrim(cellstr(num2str(currentTicks', '%d')));
%     set(ax, 'XTickLabel', labels);
%     yticks(0:10:60)
%     ytickformat('%.0f')
%     ylim([0 62])
%     set(gcf,'Name',Adr1(i))
%     set(gca, 'FontName', 'CMU Serif')
%     set(gca, "FontSize", 8)
%     if i==3
%         legend("Location","northwest")
%     else
%         legend("Location","southeast")
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified code to reduce data points to approcimately 500 for the use of
% tikz 


clear; clc; close all;

% Read the table
T = readtable("../DataComparisonStrainGaugesandForce.xlsx", Range="FA3:FI9074");

% Strain columns and number of channels per group
Adr1 = ["R0","R90","R45"; "A0","A90",""; "S2","S1",""];
Adr2 = [3, 2, 2];

% Downsample each strain column after applying threshold
for i = 1:3
    for j = 1:Adr2(i)

        % Extract strain and force
        TempX = T{:, Adr1(i,j)};
        TempY = T{:, "F"};   % Force column stays unchanged

        % Apply threshold on strain
        k = 1;
        while k <= numel(TempX) && abs(TempX(k)) < 19300
            k = k + 1;
        end
        TempX(k+1:end) = NaN;

        % Find valid indices
        valid_idx = find(~isnan(TempX));

        % Downsample indices to approx. 500 points
        N = numel(valid_idx);
        if N > 500
            idx_ds = valid_idx(round(linspace(1, N, 500)));
        else
            idx_ds = valid_idx;
        end

        % Store back downsampled strain in the table
        Xcol = nan(height(T),1);
        Xcol(idx_ds) = TempX(idx_ds);
        T{:, Adr1(i,j)} = Xcol;

    end
end

% Plotting
for i = 1:3
    figure('Color','w')
    hold on
    for j = 1:Adr2(i)
        x = T{:, Adr1(i,j)};
        y = T{:, "F"};
        valid = ~isnan(x);
        plot(x(valid), y(valid), "DisplayName", Adr1(i,j))
    end
    hold off
    xlabel("Microstrain")
    ylabel("Force [kN]")
    xtickformat('%,g')
    xticks(-20000:10000:20000)
    ax = gca; 
    currentTicks = get(ax, 'XTick');
    labels = strtrim(cellstr(num2str(currentTicks', '%d')));
    set(ax, 'XTickLabel', labels);
    yticks(0:10:60)
    ytickformat('%.0f')
    ylim([0 62])
    set(gcf,'Name',Adr1(i))
    set(gca, 'FontName', 'CMU Serif')
    set(gca, "FontSize", 8)
    if i==3
        legend("Location","northwest")
    else
        legend("Location","southeast")
    end
end
