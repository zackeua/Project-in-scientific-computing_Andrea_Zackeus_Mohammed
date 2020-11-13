clear;
clc;

% Representation av polynom i MATLAB:
% ex: a = [1,2,3] => 1*x^2 + 2*x^1 + 3*x^0,
% b = [5,0,-4,0] = 5*x^3 + 0*x^2 + (-4)*x^1 + 0*x^0
% alltså polynomen representeras av koefficienterna framför x termerna med
% högst grad först och om någon grad inte är med så är den koefficienten en
% 0:a.
%
% För att multiplicera två polynom med varandra så använder man conv() som
% tar två polynom som input och returnar ett polynom som output
%
% ex: a = [1,2,3], b = [5,0,-4,0]
% p = conv(a,b) alltså p = conv([1,2,3],[5,0,-4,0])
% p = [5, 10, 11, -8, -12, 0]
% om vi utför a*b för hand så får vi
% 1*5   1*0   1*(-4)   1*0
%       2*5   2*0      2*(-4)   2*0
%             3*5      3*0      3*(-4)   3*0
% 5    0   (-4)    0
%     10     0   (-8)     0
%           15     0   (-12)   0
% 5   10    11    -8    -12    0
% vilket vi får ut från conv(a,b)
%
% För att derivera ett polynom i MATLAB använder vi q = polyder(p)
% där p och q är två polynom
% ex: p = [1,2,3] alltså 1*x^2 + 2*x^1 + 3*x^0
% derivera för hand ger: 2*1*x^1 + 2*x^0 + 0, matlab notation: [2, 2]
% om vi kallar på q = polyder([1,2,3]) så får vi q = [2, 2]
%
% För att skriva integralen för ett polynom får man kalla på q = polyint(p)
% där p är ett polynom och q är integralen över p OBS! obestämda gränser
% för integralen nu!
%
% För att beräkna integralen mellan två punkter får man ta och evaluera
% integralen för två punkter, använd alltså values = polyval(q,[lowerBound,upperBound])
% då kommer polyval evaluera q vid lowerBound och q vid upperBound och
% returnera values = [q(lowerBound), q(upperBound)] <- tänk f(x) fast q(x)
%
% För att slutligen få värdet på integralen tar man skillnaden på
% q(upperBound)-q(lowerBound) genom diff(values)
% i koden nedan slås dessa två sista steg ihop till 1 steg med, då
% slipper vi lagra values temporärt:
% ex från koden: diff(polyval(l1,[a,b]));
% OBS! MATLABS polyint, polyder, polyval vill ha polynomen som rad matriser
% alltså a = [1, 2, 3] istället för
% a = [1;2;3] = [1;
%                2;
%                3]
% därav den nya transposen i coeff funktionen!


degree = 1; % fungerar BARA för fallet när degree = 1 (alltså P1 element) 
            % vi behöver kolla hur matriserna ska assemblas för högre grad
            % på polynomet ex. P2, P3, m.m.
h = 1;
phi = coeff(degree,h); % ta reda på koefficienterna för basfunktionerna
phiPrim = zeros(degree+1,degree);
for i=1:degree+1 % ta reda på derivatan till respektive basfunktion
    phiPrim(i,:) = polyder(phi(i,:));
end


m1 = polyint(conv(phi(1,:),phi(1,:))); % 1, 1 => j == i, skapar integralen för phi_1 * phi_1
m2 = polyint(conv(phi(2,:),phi(2,:))); % 2, 2 => j == i, likadant för de här fast andra index på phi
m3 = polyint(conv(phi(1,:),phi(2,:))); % 1, 2 => j == i±1

l1 = polyint(conv(phiPrim(1,:),phi(1,:))); % 1, 1 => j == i skapar integraler för phiPrim_1 * phi_1
l2 = polyint(conv(phiPrim(1,:),phi(2,:))); % 1, 2 => j == i-1, andra index här också men samma resonemnag
l3 = polyint(conv(phiPrim(2,:),phi(1,:))); % 2, 1 => j == i+1
l4 = polyint(conv(phiPrim(2,:),phi(2,:))); % 2, 2 => j == i

a = 0;
b = degree*h; % sätter övre gränsen till degree punkter frammåt där varje punkt är avstånd h från varandra
M1 = diff(polyval(m1,[a,b])); % beräkna integralen m1 från a till b
M2 = diff(polyval(m2,[a,b])); % liknande här
M3 = diff(polyval(m3,[a,b]));

L1 = diff(polyval(l1,[a,b])); % liknande här fast för l1 integraler o.s.v.
L2 = diff(polyval(l2,[a,b]));
L3 = diff(polyval(l3,[a,b]));
L4 = diff(polyval(l4,[a,b]));



disp(['M diagonal elements: ' num2str(M1+M2)]) % summan av de två integralerna
disp(['M lower diagonal elements: ' num2str(M3)]) % samma för upper och lower diaginal i det här fallet
disp(['M upper diagonal elements: ' num2str(M3)]) % samma för upper och lower diaginal i det här fallet

disp(['L diagonal elements: ' num2str(L1+L4)]) % summan av de två integralerna
disp(['L lower diagonal elements: ' num2str(L2)]) % när j = i-1
disp(['L upper diagonal elements: ' num2str(L3)]) % när j = i+1