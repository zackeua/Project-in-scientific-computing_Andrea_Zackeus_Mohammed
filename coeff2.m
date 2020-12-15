function c = coeff2(degree,x)
% FEM    
%{
    Function to calculate coefficients to
    lagrange basis functions for a
    given degree and step length, h.
    %}
    A = zeros(degree+1);% matris of coeff, = antalet punkter=degree+1
    for i = 1:degree+1 % varje kolonn (i) "bakifr??n" i matris A motsvarar x^degree
        A(:,i) = x.^(degree-i+1);
    end
    c = (A\eye(degree+1))'; % l??s ut koefficienterna
    % OBS!!! Nu v??ljer jag att transposa coefficient matrisen efter??t f??r att l??ttare anv??nda den med matlabs polyint, polydiff och polyval
end