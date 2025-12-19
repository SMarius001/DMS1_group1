clear; clc; close all;

T = readtable("../DataComparisonStrainGaugesandForce.xlsx",Range="CM6:EC3057");

Adr1 = ["R0";"R90";"R45";"A0";"A90";"S2";"S1"];
Adr2 = ["Test1";"Test2"];
Adr3 = ["Cycle1";"Cycle2";"Cycle3"];

for i = 1:7
    figure('Color','w')
    hold on
    for j = 1:2
        for k = 1:3
            Tempx = T.Var1;
            Tempy = T{:,string(append(Adr1(i),Adr2(j),Adr3(k)))};
            Tempy = Tempy(1:floor(length(Tempy)/200):end);
            Tempx = Tempx(1:floor(length(Tempx)/200):end);
            plot(Tempx,Tempy)
        end
    end
    hold off
    legend(Adr2(1)+Adr3(1),Adr2(1)+Adr3(2),Adr2(1)+Adr3(3),Adr2(2)+Adr3(1),Adr2(2)+Adr3(2),Adr2(2)+Adr3(3),'Location','southeast')
    xlabel("Time [s]",'Interpreter','latex')
    ylabel("Microstrain")
    xticks(0:1000:3000)
    xlim([0 3050])
    set(gcf,'Name',Adr1(i))
    set(gca, 'FontName', 'CMU Serif')
    set(gca, "FontSize", 8)
end