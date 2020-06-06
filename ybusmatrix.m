% This file contains all the impedances which are converted to admittances
% The admittance matrix is then created using those values

%generator impedances
zgen1 = 0.20*i;
zgen2 = 0.15*i;
zgen3 = 0.25*i;

%line impedances
zline14 = 0.035 + 0.225*i;
zline15 = 0.025 + 0.105*i;
zline16 = 0.040 + 0.215*i;
zline46 = 0.028 + 0.125*i;
zline56 = 0.026 + 0.175*i;

%transformer impedances
zt24 = 0.035*i;
zt35 = 0.042*i;

%generator admittances
agen1 = 1/zgen1;
agen2 = 1/zgen2;
agen3 = 1/zgen3;

%line admittances
aline14 = 1/zline14;
aline15 = 1/zline15;
aline16 = 1/zline16;
aline46 = 1/zline46;
aline56 = 1/zline56;

%transformer admittances
at24 = 1/zt24;
at35 = 1/zt35;

%--------------------------------------------------------------------------------------------------------------------------------------------------

%row 1 admittance matrix
y11 = aline14 + aline15 + aline16 + agen1;
y12 = 0;
y13 = 0;
y14 = -aline14;
y15 = -aline15;
y16 = -aline16;

%row 2 admittance matrix
y21 = 0;
y22 = agen2 + at24;
y23 = 0;
y24 = -at24;
y25 = 0;
y26 = 0;

%row 3 admittance
y31 = 0;
y32 = 0;
y33 = agen3 + at35;
y34 = 0;
y35 = -at35;
y36 = 0;

%row 4 admittance
y41 = -aline14;
y42 = -at24;
y43 = 0; 
y44 = aline14 + at24 + aline46;
y45 = 0; 
y46 = -aline46;

%row 5 admittance
y51 = -aline15;
y52 = 0;
y53 = -at35;
y54 = 0; 
y55 = aline15 + at35 + aline56;
y56 = -aline56;

%row 6 admittance
y61 = -aline16;
y62 = 0;
y63 = 0;
y64 = -aline46;
y65 = -aline56;
y66 = aline16 + aline46 + aline56;

%final admittance matrix
ybus = [ y11, y21, y13, y14, y15, y16;
y21, y22, y23, y24, y25, y26;
y31, y32, y33, y34, y35, y36;
y41, y42, y43, y44, y45, y46;
y51, y52, y53, y54, y55, y56;
y61, y62, y63, y64, y65, y66];
