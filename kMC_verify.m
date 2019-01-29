close all;
clear all;
clc;
[a,b,c,d,e,f]=ndgrid(7:9);
A=[a(:),b(:),c(:),d(:),e(:),f(:)];
whos A
for i=1:1:729;
    for j=2:1:5;
        if A(i,j)==A(i,j+1) || A(i,j)==A(i,j-1);
            A(i,:)=[6 6 6 6 6 6];
        end
    end
end
j=1
for i=1:1:length(A);
    if A(i,1)~=6;
        B(j,:)=A(i,:);
        j=j+1;
    end
end