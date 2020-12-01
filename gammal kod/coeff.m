function c = coeff(degree,h)
% FEM    
%{
    Function to calculate coefficients to
    lagrange basis functions for a
    given degree and step length, h.
    %}

    %h = 0.1;
    %degree = 1;
    % Grundtanken bakom A*C=I kommer fr??n hur vi ber??knade koefficienterna
    % C coefficientmatrisen
    % till baskunktionerna i 2d fallet i till??mpad fem kursen d??r hade vi
    % fallet
    % [1, x1, y1;  [a   [1
    %  1, x2, y2; * b =  0
    %  1, x3, y3;]  c]   0]
    % f??r den hattfunktionen som ??r 1 o punkt (x1,y1) och 0 i punkterna
    % (x2,y2) och (x3,y3), d?? kommer a vara offset i hjdled f??r
    % hattfunktionen, b ??r lutnigen i x led och c ??r lutning i y led.
    % men i v??rat 1d fall m??blerar vi om i matrisen n??got...
    % t.ex. fallet f??r grad 2 p?? hattfunktioner
    % [x1^2, x1^1, x1^0;  [a   [1
    %  x2^2, x2^1, x2^0; * b =  0
    %  x3^2, x3^1, x3^0]   c]   0]
    % h??r ??r fallet att hattfunktionen kommer vara 1 vid punkten x=x1 och 0
    % vid punkterna x=x2 och x=x3. koefficienten a ??r nu ist??llet faktorn
    % framf??r x^2 termen, b ??r framf??r x^1 allts?? x och c ??r framf??r x^0
    % allts?? 1
    % f??r att l??sa ut alla basfunktioners koefficienter direkt s?? s??tter vi
    % h??gerledet till en identitetsmatris ist??llet, vi f??r ??ven ut??ka
    % koefficient matrisen till tre kolonner f??r att f?? med de
    % funktionernas koefficienter, egentligen ??r det ingen skillnad med att
    % g??ra samma sak som ovan fast med olika h??gerled, blir ??got enklare
    % kod i MATLAB s??h??r bara.
    % [x1^2, x1^1, x1^0;  [a1, a2, a3   [1, 0, 0;
    %  x2^2, x2^1, x2^0; * b1, b2, b3 =  0, 1, 0;
    %  x3^2, x3^1, x3^0]   c1, c2, c3]   0, 0, 1]
    
    X = (0:degree) * h; % skapa en vektor med degree+1 antal x positioner fr??n 0 och fram??t med h i avst??nd fr??n varandra
                        % p?? dessa positioner vill vi l??sa polynomen till
                        % antingen 1 eller 0 (fixeras i dentitetsmatrisen
    A = zeros(degree+1);% matris of coeff, = antalet punkter=degree+1

    for i = 1:degree+1 % varje kolonn (i) "bakifr??n" i matris A motsvarar x^degree
        A(:,i) = X.^(degree-i+1);
    end
    c = (A\eye(degree+1))'; % l??s ut koefficienterna
    % OBS!!! Nu v??ljer jag att transposa coefficient matrisen efter??t f??r att l??ttare anv??nda den med matlabs polyint, polydiff och polyval
end