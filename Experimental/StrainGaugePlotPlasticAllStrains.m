clear; clc; close all;

T = readtable("../DataComparisonStrainGaugesandForce.xlsx",Range="EQ3:EX13203");

Adr1 = ["R0","R90","R45";"A0","A90","";"S2","S1",""];
Adr2 = [3, 2, 2];

for i = 1:3
    figure('Color','w')
    hold on
    for j = 1:Adr2(i)
        Tempx = T.t;
        Tempy = T{:,Adr1(i,j)};
        Tempy = Tempy(1:floor(length(Tempy)/200):end);
        Tempx = Tempx(1:floor(length(Tempx)/200):end);
        plot(Tempx,Tempy)
    end
    hold off
    xlabel("$t$ [s]",'Interpreter','latex')
    ylabel("Microstrain")
    xticks(0:400:1400)
    xlim([0 1400])
    ylim([-20000 20000])
    set(gcf,'Name',Adr1(i))
    set(gca, 'FontName', 'CMU Serif')
    set(gca, "FontSize", 8)
    legend("Location","southwest")
end