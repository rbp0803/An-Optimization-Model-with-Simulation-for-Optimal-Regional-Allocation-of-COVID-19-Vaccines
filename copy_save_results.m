function [results, sum_tanan]=copy_save_results(a,b,c,d,e,f,x,y,z)

Y1=a; Y2=b; Y3=c; Z1=d; Z2=e; Z3=f; ZZ1=x; ZZ2=y; ZZ3=z;

%primary dose
YY11=[];
YY12=[];
YY13=[];
YY14=[];
YY15=[];

   for ik=1:17
       YY11=cat(1,YY11,Y1(5*(ik-1)+1,:));
       YY12=cat(1,YY12,Y1(5*(ik-1)+2,:));
       YY13=cat(1,YY13,Y1(5*(ik-1)+3,:));
       YY14=cat(1,YY14,Y1(5*(ik-1)+4,:));
       YY15=cat(1,YY15,Y1(5*ik,:));
   end
   
 YY21=[];
 YY22=[];
 YY23=[];
       for ik=1:17
       YY21=cat(1,YY21,Y2(3*(ik-1)+1,:));
       YY22=cat(1,YY22,Y2(3*(ik-1)+2,:));
       YY23=cat(1,YY23,Y2(3*ik,:));
       end
       
 YY31=[];
 YY32=[];
 YY33=[];
       for ik=1:17
       YY31=cat(1,YY31,Y3(3*(ik-1)+1,:));
       YY32=cat(1,YY32,Y3(3*(ik-1)+2,:));
       YY33=cat(1,YY33,Y3(3*ik,:));
       end
       
%booster       
ZZ11=[];
ZZ12=[];
ZZ13=[];
ZZ14=[];
ZZ15=[];

   for ik=1:17
       ZZ11=cat(1,ZZ11,Z1(5*(ik-1)+1,:));
       ZZ12=cat(1,ZZ12,Z1(5*(ik-1)+2,:));
       ZZ13=cat(1,ZZ13,Z1(5*(ik-1)+3,:));
       ZZ14=cat(1,ZZ14,Z1(5*(ik-1)+4,:));
       ZZ15=cat(1,ZZ15,Z1(5*ik,:));
   end       
 

 ZZ21=[];
 ZZ22=[];
 ZZ23=[];
       for ik=1:17
       ZZ21=cat(1,ZZ21,Z2(3*(ik-1)+1,:));
       ZZ22=cat(1,ZZ22,Z2(3*(ik-1)+2,:));
       ZZ23=cat(1,ZZ23,Z2(3*ik,:));
       end   
   
 ZZ31=[];
 ZZ32=[];
 ZZ33=[];
       for ik=1:17
       ZZ31=cat(1,ZZ31,Z3(3*(ik-1)+1,:));
       ZZ32=cat(1,ZZ32,Z3(3*(ik-1)+2,:));
       ZZ33=cat(1,ZZ33,Z3(3*ik,:));
       end
       
%second booster       
Z2Z11=[];
Z2Z12=[];
Z2Z13=[];
Z2Z14=[];
Z2Z15=[];

   for ik=1:17
       Z2Z11=cat(1,Z2Z11,ZZ1(5*(ik-1)+1,:));
       Z2Z12=cat(1,Z2Z12,ZZ1(5*(ik-1)+2,:));
       Z2Z13=cat(1,Z2Z13,ZZ1(5*(ik-1)+3,:));
       Z2Z14=cat(1,Z2Z14,ZZ1(5*(ik-1)+4,:));
       Z2Z15=cat(1,Z2Z15,ZZ1(5*ik,:));
   end

 Z2Z21=[];
 Z2Z22=[];
 Z2Z23=[];
       for ik=1:17
       Z2Z21=cat(1,Z2Z21,ZZ2(3*(ik-1)+1,:));
       Z2Z22=cat(1,Z2Z22,ZZ2(3*(ik-1)+2,:));
       Z2Z23=cat(1,Z2Z23,ZZ2(3*ik,:));
       end

Z2Z31=[];
Z2Z32=[];
Z2Z33=[];
       for ik=1:17
       Z2Z31=cat(1,Z2Z31,ZZ3(3*(ik-1)+1,:));
       Z2Z32=cat(1,Z2Z32,ZZ3(3*(ik-1)+2,:));
       Z2Z33=cat(1,Z2Z33,ZZ3(3*ik,:));
       end       
       
   
       
results=cell(11, 3);
results{1,1}=YY11; results{1,2}=ZZ11; results{1,3}=Z2Z11;
results{2,1}=YY12; results{2,2}=ZZ12; results{2,3}=Z2Z12;
results{3,1}=YY13; results{3,2}=ZZ13; results{3,3}=Z2Z13;
results{4,1}=YY14; results{4,2}=ZZ14; results{4,3}=Z2Z14;
results{5,1}=YY15; results{5,2}=ZZ15; results{5,3}=Z2Z15;
results{6,1}=YY21; results{6,2}=ZZ21; results{6,3}=Z2Z21;
results{7,1}=YY22; results{7,2}=ZZ22; results{7,3}=Z2Z22;
results{8,1}=YY23; results{8,2}=ZZ23; results{8,3}=Z2Z23;       
results{9,1}=YY31; results{9,2}=ZZ31; results{9,3}=Z2Z31;
results{10,1}=YY32; results{10,2}=ZZ32; results{10,3}=Z2Z32;
results{11,1}=YY33; results{11,2}=ZZ33; results{11,3}=Z2Z33;

     
sum_tanan=[sum(sum([YY11 YY12 YY13 YY14 YY15 ZZ11 ZZ12 ZZ13 ZZ14 ZZ15 Z2Z11 Z2Z12 Z2Z13 Z2Z14 Z2Z15...
    YY21 YY22 YY23 ZZ21 ZZ22 ZZ23 Z2Z21 Z2Z22 Z2Z23 YY31 YY32 YY33 ZZ31 ZZ32 ZZ33 Z2Z31 Z2Z32 Z2Z33])); ...
    sum(sum([YY11 ZZ11 Z2Z11 YY21 ZZ21 Z2Z21 YY31 ZZ31 Z2Z31])); sum(sum([YY12 ZZ12 Z2Z12 YY22 ZZ22 Z2Z22  YY32 ZZ32 Z2Z32]));...
    sum(sum([YY13 ZZ13 Z2Z13])); sum(sum([YY14 ZZ14 Z2Z14]));...
    sum(sum([YY15 ZZ15 Z2Z15 YY23 ZZ23 Z2Z23 YY33 ZZ33 Z2Z33]))];


