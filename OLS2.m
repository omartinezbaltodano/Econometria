function [b R2] = OLS2(X,y,c)
%-----------------------------------------------
% PROPOSITO: Hace la estimacion de los parametros por OLS
%-----------------------------------------------
% INSUMOS  : X  : nxK matriz de variables independientes (sin constante)
%            y  : nx1 vector de variable dependiente
%            c  : 1x1 : 1 si hay constante en el modelo
%                       0 si no hay constante en el modelo
%-----------------------------------------------
% OUTPUT   : b  : Kx1 vector de parametros 
%            R2 : 1x1 coeficiente de determinacion
%-----------------------------------------------
if nargin < 3  %% Asumo que en el X no hay constante y la agrego
    X = [ones(size(X,1),1) X];
    b = (X'*X)^(-1)*X'*y;
    %coeficiente de determinacion 
    n = size(y,1);
    Aux1 = X*b-((1/n)*sum(y));
    Aux2 = Aux1'*Aux1;
    Aux3 = y-((1/n)*sum(y));
    Aux4 = Aux3'*Aux3;
    R2 = Aux2/Aux4;
else 
    if c == 1
        X = [ones(size(X,1),1) X];
        b = (X'*X)^(-1)*X'*y;
        %coeficiente de determinacion 
        n = size(y,1);
        Aux1 = X*b-((1/n)*sum(y));
        Aux2 = Aux1'*Aux1;
        Aux3 = y-((1/n)*sum(y));
        Aux4 = Aux3'*Aux3;
        R2 = Aux2/Aux4;
    elseif c == 0
        b = (X'*X)^(-1)*X'*y;
        %coeficiente de determinacion no centrado
        n    = size(y,1);
        yhat = X*b;
        R2 = (yhat'*yhat)/(y'*y);
    else 
        error('c solo puedo tomar los valores {0,1}')
    end 
end 

end 