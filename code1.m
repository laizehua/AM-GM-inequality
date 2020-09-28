%{
This is the supplement code for "Recht¨CR?e Noncommutative
Arithmetic-Geometric Mean Conjecture is False."

This program requires SeDuMi in Matlab. Download SeDuMi from http://sedumi.ie.lehigh.edu/?page_id=58
Add the SeDumi folder into Matlab path.
%}

%{
VarNum is the number of matrix, n.
m is the degree of our inequality.
postive = 1 or -1. Minimize \lambda_1 and \lambda_2 respectively.
%}
tic

VarNum = 5;
m = 5;
Deg = floor(m/2);
postive = 1;

Dim1 = NumVar(Deg*2+2,VarNum);
Dim2 = NumVar(Deg*2+1,VarNum);
Dim3 = NumVar(Deg+1,VarNum);

%{
Prod is the product table for our noncommutative monomials.
%}
Prod = zeros(Dim2, Dim2);
for i = 1:Dim2
    Prod(1,i)=i;
    Prod(i,1)=i;
end

for i = 1:(Deg*2)
    for j = 1:(Deg*2)
        for x = 1:VarNum^i
            for y = 1:VarNum^j
                Prod(NumVar(i,VarNum)+x,NumVar(j,VarNum)+y) = NumVar(i+j,VarNum)+(x-1)*VarNum^j+y;
            end
        end
    end
end
%{
Tra is the transport operation of our noncommutative monomials.
%}
Tra = zeros(1,Dim3);
for i = 1:VarNum+1
    Tra(i)=i;
end
for i = 1:VarNum^2
    p = ceil(i/VarNum);
    q = i-(p-1)*VarNum;
    Tra(VarNum+1+i)=(q-1)*VarNum+p+VarNum+1;
end


%{
Const and Equ is the linear constraints. Const*x = Equ.
%}
Const = zeros(Dim1, (VarNum+1)*Dim3*Dim3+1);

for i = 0:Deg
    for j = 0:Deg
        for x = 1:VarNum^i
            for y = 1:VarNum^j
                a = Prod(Tra(NumVar(i,VarNum)+x),NumVar(j,VarNum)+y);
                b = (NumVar(i,VarNum)+x-1)*Dim3+NumVar(j,VarNum)+y;
                for k = 1:VarNum
                    a = Prod(Prod(Tra(NumVar(i,VarNum)+x),k+1),NumVar(j,VarNum)+y);
                    Const(a,b+(k-1)*Dim3^2)=1;
                    Const(a,b+VarNum*Dim3^2)=-1;
                end
                a = Prod(Tra(NumVar(i,VarNum)+x),NumVar(j,VarNum)+y);
                Const(a,b+VarNum*Dim3^2)=VarNum;
            end
        end
    end
end
Const(1, (VarNum+1)*Dim3*Dim3+1)=-1;

Equ = zeros(Dim1,1);
A = perms(1:VarNum);
A = unique(A,'rows');
for i=1:length(A)
    index = NumVar(m,VarNum);
    for j = 1:m
        index = index + (A(i,j)-1)*VarNum^(m-j);
    end
    Equ(index+1) = -postive;
end

%{
K is the dimension of SDP variables.
%}
K.s = zeros(VarNum+2,1);
for i=1:VarNum+1
    K.s(i) = Dim3;
end
K.s(VarNum+2) = 1;

%{
Tr is the variable we are minimize: \lambda
%}
Tr = zeros((VarNum+1)*Dim3*Dim3+1,1);
Tr((VarNum+1)*Dim3*Dim3+1)=1;
pars.eps=0;

%{
Call SeDuMi.

X(1:Dim3^2), ...
,X((VarNum*Dim3^2:(VarNum+1)*Dim3^2) are the sum of square matrices.

X((VarNum+1)*Dim3*Dim3+1) is \lambda
%}
[X,y,info] = sedumi(Const,Equ,Tr,K,pars);

toc