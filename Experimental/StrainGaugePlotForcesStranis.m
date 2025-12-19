clear; clc; close all;

T = readtable("../DataComparisonStrainGaugesandForce.xlsx",Range="CM6:EI3057");

Adr1 = ["R0","R90","R45";"A0","A90","";"S2","S1",""];
Adr2 = ["Test1";"Test2"];
Adr3 = ["Cycle1";"Cycle2";"Cycle3"];
var1 = [3, 2, 2];
Var2 = ["o","x"];
Var3 = [[0.0660,0.4430,0.7450];[0.8660,0.3290,0];[0.9290,0.6940,0.1250]];

for i = 1:3
    figure('Color','w')
    hold on
    h(1) = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0], 'MarkerSize', 8);
    h(2) = plot(NaN, NaN, 'x', 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0], 'MarkerSize', 8);
    for j = 1:2
        Tempx = zeros(3051,3);
        for l = 1:var1(i)
            for k = 1:3
                TempX = T{:,string(append('Force1',Adr2(j),Adr3(k)))};
                TempY = T{:,string(append(Adr1(i,l),Adr2(j),Adr3(k)))};
                TempY = TempY(1:floor(length(TempY)/200):end);
                TempX = TempX(1:floor(length(TempX)/200):end);
                s = scatter(TempY,TempX,Marker=Var2(j));
                s.MarkerFaceColor = Var3(k,:);
                s.MarkerEdgeColor = Var3(k,:);
                s.MarkerFaceAlpha = 0.8;
                s.MarkerEdgeAlpha = 0.8;
                s.LineWidth = 0.5;
                h(k+2) = plot(NaN, NaN, 'o', 'MarkerFaceColor', Var3(k,:), 'MarkerEdgeColor', Var3(k,:), 'MarkerSize', 8);
            end
            Tempx(:,l) = T{:,append(Adr1(i,l),Adr2(j),Adr3(k))};
        end
    end
    Tempy = 3;
    TempxMax = max(abs(Tempx),[],"all");
    for l = 1:var1(i)
        if i == 3
            text((Tempx(floor(end/6),l))*1.5+50,Tempy,Adr1(i,l),FontSize=8);
        else
            text((Tempx(floor(end/6),l))*2-10,Tempy,Adr1(i,l),FontSize=8);
        end
    end
    hold off
    legend(h,{'Test 1', 'Test 2', 'Cycle 1', 'Cycle 2', 'Cycle 3'},'Location', 'southeast',Box='off')
    ylabel("Force [kN]")
    xlabel("Microstrain")
    yticks(0:0.5:4)
    ytickformat('%.1f')
    %ylim([0 4.1])
    set(gcf,'Name',Adr1(i))
    set(gca, 'FontName', 'CMU Serif')
    set(gca, "FontSize", 8)
end