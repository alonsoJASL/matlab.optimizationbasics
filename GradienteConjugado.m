function [x,k]=GradienteConjugado(A,b, Tol)
%                   GRADIENTE CONJUGADO
% Encuentra la solucion para el sistema de ecuaciones lineales
%                        Ax = b
% por el metodo de gradiente conjugado. 
% 
% usage:     [x,k]=GradienteConjugado(A,b)
%         donde: 
%        SALIDA
%           x := solucion al sistema de ecuaciones
%           k := numero de iteraciones que toma el metodo
%        ENTRADA
%           A := Matriz simetrica y positiva definida.
%           b := vector de soluciones al sistema.
% 
% 

%
% Jose Alonso Solis Lemus (2012. ITAM) 
% for license and more code check: 
% https://github.com/alonsoJASL/matlab.optimizationbasics.git
%

n=length(b);
x=zeros(n,1);
if nargin < 3
    tol=1.e-06;
else
    tol=Tol;
end

rk=-b; %valor de r0
pk=-rk;   %valor de p0

k=0; %el contador se inicia en cero

while (norm(rk)>tol)
    den=pk'*A*pk;
    alfak=-(rk'*pk)/(den);
    x=x+alfak*pk; %calculo nuevo punto x
    rk=A*x-b;% con el punto de la linea anterior calculo el nuevo residual
    betak= (rk'*A*pk)/ (den); 
    pk=-rk+betak*pk;
    
    k=k+1;  
end

end