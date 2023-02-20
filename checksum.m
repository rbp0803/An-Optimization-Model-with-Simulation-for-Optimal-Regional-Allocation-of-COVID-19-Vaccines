function [a_05, a_12, a_18, a_05_boost, a_12_boost, a_18_boost, a_05_boost2, a_12_boost2, a_18_boost2]=checksum(a)
%per week 
%a=result

%primary
a18=a(1:5,1);
a12=a(6:8,1);
a05=a(9:11,1);

%18 onwards
aa3=[];
for i=1:length(a18)
   if i==3
   aa=a18(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
       
   else
   aa=a18(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa2=aa2*0.5;
   aa3=cat(2,aa3,aa2);
   
   end
end
aa3=aa3; 

a_18=sum(aa3,2);

%12to17
aa3=[];
for i=1:length(a12)
   aa=a12(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
end
   aa3=aa3*0.5;

a_12=sum(aa3,2);

%05to11
aa3=[];
for i=1:length(a05)
   aa=a05(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);

end
   aa3=aa3*0.5;
a_05=sum(aa3,2);

%first booster
a18=a(1:5,2);
a12=a(6:8,2);
a05=a(9:11,2);

aa3=[];
for i=1:length(a18)
   aa=a18(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);

end

a_18_boost=sum(aa3,2);


%12to17
aa3=[];
for i=1:length(a12)
   aa=a12(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
end

a_12_boost=sum(aa3,2);

%05to11
aa3=[];
for i=1:length(a05)
   aa=a05(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
end
a_05_boost=sum(aa3,2);


%second booster
a18=a(1:5,3);
a12=a(6:8,3);
a05=a(9:11,3);

aa3=[];
for i=1:length(a18)
   aa=a18(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);

end

a_18_boost2=sum(aa3,2);


%12to17
aa3=[];
for i=1:length(a12)
   aa=a12(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
end

a_12_boost2=sum(aa3,2);

%05to11
aa3=[];
for i=1:length(a05)
   aa=a05(i);
   aa1=cell2mat(aa);
   aa2=sum(aa1,2);
   aa3=cat(2,aa3,aa2);
end
a_05_boost2=sum(aa3,2);

end