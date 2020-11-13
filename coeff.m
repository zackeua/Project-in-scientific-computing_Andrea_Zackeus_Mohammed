function c = coeff(degree,h)
    %{
    Function to calculate coefficients to
    lagrange basis functions for a
    given degree and step length, h.
    %}

    %h = 0.1;
    %degree = 1;
    % Grundtanken bakom A*C=I kommer från hur vi beräknade koefficienterna
    % till baskunktionerna i 2d fallet i tillämpad fem kursen där hade vi
    % fallet
    % [1, x1, y1;  [a   [1
    %  1, x2, y2; * b =  0
    %  1, x3, y3;]  c]   0]
    % för den hattfunktionen som är 1 o punkt (x1,y1) och 0 i punkterna
    % (x2,y2) och (x3,y3), då kommer a vara offset i hjdled för
    % hattfunktionen, b är lutnigen i x led och c är lutning i y led.
    % men i vårat 1d fall möblerar vi om i matrisen något...
    % t.ex. fallet för grad 2 på hattfunktioner
    % [x1^2, x1^1, x1^0;  [a   [1
    %  x2^2, x2^1, x2^0; * b =  0
    %  x3^2, x3^1, x3^0]   c]   0]
    % här är fallet att hattfunktionen kommer vara 1 vid punkten x=x1 och 0
    % vid punkterna x=x2 och x=x3. koefficienten a är nu istället faktorn
    % framför x^2 termen, b är framför x^1 alltså x och c är framför x^0
    % alltså 1
    % för att lösa ut alla basfunktioners koefficienter direkt så sätter vi
    % högerledet till en identitetsmatris istället, vi får även utöka
    % koefficient matrisen till tre kolonner för att få med de
    % funktionernas koefficienter, egentligen är det ingen skillnad med att
    % göra samma sak som ovan fast med olika högerled, blir ågot enklare
    % kod i MATLAB såhär bara.
    % [x1^2, x1^1, x1^0;  [a1, a2, a3   [1, 0, 0;
    %  x2^2, x2^1, x2^0; * b1, b2, b3 =  0, 1, 0;
    %  x3^2, x3^1, x3^0]   c1, c2, c3]   0, 0, 1]
    
    X = (0:degree) * h; % skapa en vektor med degree+1 antal x positioner från 0 och framåt med h i avstånd från varandra
                        % på dessa positioner vill vi låsa polynomen till
                        % antingen 1 eller 0 (fixeras i dentitetsmatrisen
    A = zeros(degree+1);

    for i = 1:degree+1 % varje kolonn (i) "bakifrån" i matris A motsvarar x^degree
        A(:,i) = X.^(degree-i+1);
    end
    c = (A\eye(degree+1))'; % lös ut koefficienterna
    % OBS!!! Nu väljer jag att transposa coefficient matrisen efteråt för att lättare använda den med matlabs polyint, polydiff och polyval
end