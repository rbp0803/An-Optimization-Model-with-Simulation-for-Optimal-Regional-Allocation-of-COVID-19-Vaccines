function [results, sum_tanan]=copy_save_results(a,b,c,d,e,f)

Y1=a; Y2=b; Y3=c; Z1=d; Z2=e; Z3=f;

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
       for ik=1:17
       YY31=cat(1,YY31,Y3(2*(ik-1)+1,:));
       YY32=cat(1,YY32,Y3(2*ik,:));
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
       for ik=1:17
       ZZ31=cat(1,ZZ31,Z3(2*(ik-1)+1,:));
       ZZ32=cat(1,ZZ32,Z3(2*ik,:));
       end       
       
results=cell(10, 2);
results{1,1}=YY11; results{1,2}=ZZ11;
results{2,1}=YY12; results{2,2}=ZZ12;
results{3,1}=YY13; results{3,2}=ZZ13;
results{4,1}=YY14; results{4,2}=ZZ14;
results{5,1}=YY15; results{5,2}=ZZ15;
results{6,1}=YY21; results{6,2}=ZZ21;
results{7,1}=YY22; results{7,2}=ZZ22;
results{8,1}=YY23; results{8,2}=ZZ23;       
results{9,1}=YY31; results{9,2}=ZZ31;
results{10,1}=YY32; results{10,2}=ZZ32;

     
sum_tanan=[sum(sum([YY11 YY12 YY13 YY14 YY15 ZZ11 ZZ12 ZZ13 ZZ14 ZZ15 YY21 YY22 YY23 ZZ21 ZZ22 ZZ23 YY31 YY32 ZZ31 ZZ32])); ...
    sum(sum([YY11 ZZ11 YY21 ZZ21 YY31 ZZ31])); sum(sum([YY12 ZZ12 YY22 ZZ22]));...
    sum(sum([YY13 ZZ13])); sum(sum([YY14 ZZ14]));...
    sum(sum([YY15 ZZ15 YY23 ZZ23 YY32 ZZ32]))];


