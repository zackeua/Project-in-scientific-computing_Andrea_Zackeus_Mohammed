function w = weights(p,m)
% p: endast j?mna pga central FDM
    %{
    Function to generate weights, w for a derivative of order m,
    in our case m=1, central differnence stencil for order p.
    found at: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    %}
    p = p/2; % order of central difference stencil
    %m = 1; % order of derivative
    
    % algoritmen ??r stulen direkt fr??n wikipedial??nken ovan
    A = ones(2*p+1);
    b = zeros(2*p+1,1);
    ROW = -p:p; % varje rad best??r av talen mellan -p till p upph??jt till radnumret med start p?? 0
    b(m+1) = factorial(m); % man s??tter rad (m+1) i b vektorn till m!
    for i=1:2*p+1
        A(i,:) = ROW.^(i-1); % s??tt rad nummer (i) till ROW.^(i-1), det h??r ??r ekvivalent med loopen nedanf??r
        % for j=1:2*p+1
        %     A(i,j) = A(i,j)*(-p+j-1)^(i-1);
        % end

    end
    w = (A\b)'; % vi f??r inte gl??mma att dividera med h^m efter/medan vi assemblar matrisen
end