% -------------------------------------------------------------------------------------%
% Optimal Allocation of Covid Vaccine Brands in Different Regions in the Philippines   %                           
% by Rodney Pino, Renier Mendoza, Erika Antonette Enriquez, Victoria May Mendoza, and  %
% Arriane Crystal Velasco                                                              % 
% -------------------------------------------------------------------------------------%

clear
tic;
N=17; %Number of Regions              
O=49; %Number of Weeks 

MM=[5,3,3]; %Number of authorized vaccine brands per age group, [18 and up, 12-17, 05-11].

Vsup=[5920929 171506 0 497635 3907209]';  %Available Vaccine Supplies at the start of simulation. In order - [Pfizer Moderna Janssen Astrazeneca Sinovac] 

[~, ~, Irate]=xlsread('Incidence_rate.xlsx'); % Incidence Rate - Proportion of New Cases 
Irate=Irate(1:N,1);                          % to the number of Population per region in
Irate=cell2mat(Irate);                        % the last 14 Days before Jan 20, 2023  




%% Objective Function Setup
%The objective function represents the Covid-19 infections after vaccine
%inoculation.
%===================================For 18 onwards===================================%
%Authorized Vaccines Arranged in Order
% [Pfizer, Moderna, Janssen, Astra/Sputnik V, Sinovac]
e=[.640;.440;0.470;0.498;0.269];    %Effectiveness Against Omicron 
ew=[.470;.235;.37;.304;.05];      %Waning Effect Against Omicron
eb=[.770;.721;.54;.556;0.15];     %Booster Effectiveness Against Omicron
eww=[.53;.512;.37;.286;.004];
ebb=[.770;.721;.540;.556;.15];    %Second Booster Effectiveness
M=5;
f=zeros(M*O*N,1);
fb=f; fbb=fb;
for k=1:O
    for j=1:N
    for i=1:M
        if (i==3) 
            eta=0.5; 
        else
            eta=1;
        end
            if k<=O-24
        f((k-1)*N*M+(j-1)*M+i)=(-e(i)/(2*eta))*Irate(j);
            end
        
           if k<=O-12 
        fb((k-1)*N*M+(j-1)*M+i)=(1-ew(i))*(1-eb(i))*Irate(j);
           end
           
        fbb((k-1)*N*M+(j-1)*M+i)=(1-eww(i))*(1-ebb(i))*Irate(j);
    end
    end
end
%====================================================================================%

%====================================For 12 to 17====================================%
% Authorized Vaccines Arranged in Order
%[Pfizer, Moderna, Sinovac]
e=[.370;.379;.269];         %Effectiveness Against Omicron
ew=[.266;.268;.05];       %Waning Effect Against Omicron
eb=[.479;0.477;0];        %Booster Effectiveness Against Omicron ; Only Pfizer and Moderna is authorized for booster
eww=[.41;.435;0];
ebb=[.479;0.477;0];      %Second Booster Effectiveness   
M=3;
f1=zeros(M*O*N,1);
fb1=f1; fbb1=fb1;
eta=1;
for k=1:O
    for j=1:N
    for i=1:M
                if k<=O-24
        f1((k-1)*N*M+(j-1)*M+i)=(-e(i)/(2*eta))*Irate(j);
                end
                
                if k<=O-12
        fb1((k-1)*N*M+(j-1)*M+i)=(1-ew(i))*(1-eb(i))*Irate(j);
                end
        fbb1((k-1)*N*M+(j-1)*M+i)=(1-eww(i))*(1-ebb(i))*Irate(j);
    end
    end
end
%====================================================================================%

%====================================For 05 to 11====================================%
% Authorized Vaccines Arranged in Order
%[Pfizer, Moderna, Sinovac]
e=[.49;0.48;.269];        %Effectiveness Against Omicron
ew=[0.27;0.41;.05];     %Waning Effect Against Omicron
eb=[0.58;0.61;0];         %Booster Effectiveness Against Omicron; Set to Pfizer and Moderna authorized for booster
eww=[0.24;0.33;0];
ebb=[0.58;0.61;0];        %Second Booster Effectiveness
M=3;
f2=zeros(M*O*N,1);
fb2=f2; fbb2=fb2;
eta=1;
for k=1:O
    for j=1:N
    for i=1:M
                if k<=O-24
        f2((k-1)*N*M+(j-1)*M+i)=(-e(i)/(2*eta))*Irate(j);
                end
                
                if k<=O-12
        fb2((k-1)*N*M+(j-1)*M+i)=(1-ew(i))*(1-eb(i))*Irate(j);
                end
        fbb2((k-1)*N*M+(j-1)*M+i)=(1-eww(i))*(1-ebb(i))*Irate(j);
    end
    end
end

%====================================================================================%
coeff=[f; fb; fbb; f1; fb1; fbb1; f2; fb2; fbb2]; %Objective Function Coefficients

%% Constraint on Vax supply
%Apparently, the supply of covid vaccines in the Philippines is insufficient 
%to inoculate all the eligible population for full primary vaccine series, first 
%booster jab, and second booster shot. Thus, the current supply of each vaccine brand is set to be the
%minimum for this constraint. 
%===================================For 18 onwards===================================%
A_vax_supply=[];
M=5;
for i=1:M
    C=zeros(1,M*O*N);
    Cc=C;
    Cb=C;
    Cb2=C;
    for k=1:O
    
    for j=1:N

            if k<=O-24
        C((k-1)*N*M+(j-1)*M+i)=1;
            end
            if k<=O-12
        Cb((k-1)*N*M+(j-1)*M+i)=1;
            end
        Cb2((k-1)*N*M+(j-1)*M+i)=1;
    end
    end
    A_vax_supply=cat(1,A_vax_supply, [C Cb Cb2]);
end
%====================================================================================%
A=A_vax_supply;

%====================================For 12 to 17====================================%
A_vax_supply=[];
M=3;
for i=1:M
    C=zeros(1,M*O*N);
    Cc=C;
    Cb=C;
    Cb2=C;
    for k=1:O
    
    for j=1:N
            if k<=O-24
        C((k-1)*N*M+(j-1)*M+i)=1;
            end
            
            if k<=O-12
        Cb((k-1)*N*M+(j-1)*M+i)=1;
            end
        Cb2((k-1)*N*M+(j-1)*M+i)=1;
    end
    end
    if i==3
       Cb=zeros(1,M*O*N);   %Sinovac not yet applicable for booster
       Cb2=zeros(1,M*O*N);  %Sinovac not yet applicable for second booster
    end
    
    A_vax_supply=cat(1,A_vax_supply, [C Cb Cb2]);
    
    if i==2
        A_vax_supply=cat(1,A_vax_supply, zeros(2,M*O*N*3));
    end
end
%====================================================================================%
A=cat(2,A,A_vax_supply);


%====================================For 05 to 11====================================%
A_vax_supply=[];
M=3;
for i=1:M
    C=zeros(1,M*O*N);
    Cc=C;
    Cb=C;
    Cb2=C;
    for k=1:O
    
    for j=1:N
            if k<=O-24
        C((k-1)*N*M+(j-1)*M+i)=1;
            end
            
            if k<=O-12
        Cb((k-1)*N*M+(j-1)*M+i)=1;
            end
        Cb2((k-1)*N*M+(j-1)*M+i)=1;
    end
    
    if i==3
       Cb=zeros(1,M*O*N); %Sinovac not yet allowed for booster 
       Cb2=zeros(1,M*O*N); %Sinovac not yet applicable for second booster
    end
    
    end
    A_vax_supply=cat(1,A_vax_supply, [C Cb Cb2]);
    
    if i==2
        A_vax_supply=cat(1,A_vax_supply, zeros(2,M*O*N*3));
    end
end

%====================================================================================%
A=cat(2,A,A_vax_supply);
A=-A;
b=-Vsup;


%% Constraint on Max Population Vaccination in each Region
%The goal is to inoculate all eligible population for full primary series,
%first booster shot, and second booster vaccination. With the "Constraint on Eligible for Booster," primary 
%series vaccination is set to be done on the O-24 week as the model aims
%also to give all eligible population a second booster shot (24 weeks after complete primary
%vaccination) until the end of the simulation.
%===================================For 18 onwards===================================%

M=5;
[~, ~, raw]=xlsread('Constraint_Data_Pop_Jan20.xlsx');
Pop_No_Vaxx_18=raw(1:N,12);
Pop_No_Vaxx_18=cell2mat(Pop_No_Vaxx_18);
Pop_No_VaxxB_18=raw(1:N,16);
Pop_No_VaxxB_18=cell2mat(Pop_No_VaxxB_18); Pop_No_VaxxB_18=round(Pop_No_VaxxB_18);
Pop_No_VaxxBB_18=raw(20:36,16);
Pop_No_VaxxBB_18=cell2mat(Pop_No_VaxxBB_18); Pop_No_VaxxBB_18=round(Pop_No_VaxxBB_18);

for k=1:O-24
Cc1=[];

for j=1:N
    if k==1
    Cc=zeros(1,M*O*N);
    else
    Cc=Cbb3(j,:); 
    end
    
    for i=1:M
        if (i==3) 
            eta=0.5; 
        else
            eta=1;
        end
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    
    Cc1=cat(1,Cc1,Cc);
end

Cbb3=Cc1;

CC1_18=Cc1;

A=cat(1,A,[Cc1 0*Cc1 0*Cc1 zeros(17,MM(2)*O*N*3) zeros(17,MM(3)*O*N*3)]);
b=cat(1,b,Pop_No_Vaxx_18);

end

%============Max booster doses for 18 onwards at the end of simulation============%
Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O-12
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end
AB=[0*Cc1 Cc1 0*Cc1 zeros(17,MM(2)*O*N*3) zeros(17,MM(3)*O*N*3)];
beq=Pop_No_VaxxB_18;
%====================================================================================%
AB1=[]; beq1=[];

% A=[A;AB];
% b=[b;beq];

AB1=cat(1,AB1,AB);
beq1=cat(1,beq1,beq);


%============Second booster doses for 18 onwards at the end of simulation============%
Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end
ABB=[0*Cc1 0*Cc1 Cc1 zeros(17,MM(2)*O*N*3) zeros(17,MM(3)*O*N*3)];
beqbb=Pop_No_VaxxBB_18;
%====================================================================================%
ABB1=[]; beqq1=[];

ABB1=cat(1,ABB1,ABB);
beqq1=cat(1,beqq1,beqbb);




%====================================For 12 to 17====================================%

M=3;
eta=1;
[~, ~, raw]=xlsread('Constraint_Data_Pop_Jan20.xlsx');
Pop_No_Vaxx_18=raw(1:N,11);
Pop_No_Vaxx_18=cell2mat(Pop_No_Vaxx_18);
Pop_No_VaxxB_18=raw(1:N,15);
Pop_No_VaxxB_18=cell2mat(Pop_No_VaxxB_18); Pop_No_VaxxB_18=round(Pop_No_VaxxB_18);
Pop_No_VaxxBB_12=raw(20:36,15);
Pop_No_VaxxBB_12=cell2mat(Pop_No_VaxxBB_12); Pop_No_VaxxBB_12=round(Pop_No_VaxxBB_12);

for k=1:O-24
    Cc1=[];
for j=1:N
    if k==1
    Cc=zeros(1,M*O*N);
    else
    Cc=Cbb3(j,:); 
    end
 
    for i=1:M
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    
    Cc1=cat(1,Cc1,Cc);
end

Cbb3=Cc1;

CC1_12=Cc1;

A=cat(1,A,[zeros(17,MM(1)*O*N*3) Cc1 0*Cc1 0*Cc1 zeros(17,MM(3)*O*N*3)]);
b=cat(1,b,Pop_No_Vaxx_18); 

end

%============Max booster doses for ages 12 to 17 at the end of simulation============%

Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O-12
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end

AB=[zeros(17,MM(1)*O*N*3) 0*Cc1 Cc1 0*Cc1 zeros(17,MM(3)*O*N*3)];
beq=Pop_No_VaxxB_18;
%====================================================================================%
% A=[A;AB];
% b=[b;beq];

AB1=cat(1,AB1,AB);
beq1=cat(1,beq1,beq);

%============Second booster doses for ages 12 to 17 at the end of simulation============%

Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end

ABB=[zeros(17,MM(1)*O*N*3) 0*Cc1 0*Cc1 Cc1 zeros(17,MM(3)*O*N*3)];
beqbb=Pop_No_VaxxBB_12;
%====================================================================================%
ABB1=cat(1,ABB1,ABB);
beqq1=cat(1,beqq1,beqbb);




%====================================For 05 to 11====================================%

M=3;
eta=1;
[~, ~, raw]=xlsread('Constraint_Data_Pop_Jan20.xlsx');
Pop_No_Vaxx_18=raw(1:N,10);
Pop_No_Vaxx_18=cell2mat(Pop_No_Vaxx_18);
Pop_No_VaxxB_18=raw(1:N,14);
Pop_No_VaxxB_18=cell2mat(Pop_No_VaxxB_18); Pop_No_VaxxB_18=round(Pop_No_VaxxB_18);
Pop_No_VaxxBB_05=raw(20:36,14);
Pop_No_VaxxBB_05=cell2mat(Pop_No_VaxxBB_05); Pop_No_VaxxBB_05=round(Pop_No_VaxxBB_05);


for k=1:O-24
    Cc1=[];
for j=1:N
    if k==1
    Cc=zeros(1,M*O*N);
    else
    Cc=Cbb3(j,:);
    end
 
    for i=1:M
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
 
    Cc1=cat(1,Cc1,Cc);
end

Cbb3=Cc1;

CC1_05=Cc1;

A=cat(1,A,[zeros(17,MM(1)*O*N*3) zeros(17,MM(2)*O*N*3) Cc1 0*Cc1 0*Cc1]);
b=cat(1,b,Pop_No_Vaxx_18);

end

%============Max booster doses for ages 05 to 11 at the end of simulation============%
Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O-12
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end


AB=[zeros(17,MM(1)*O*N*3) zeros(17,MM(2)*O*N*3) 0*Cc1 Cc1 0*Cc1];
beq=Pop_No_VaxxB_18;
%====================================================================================%
% A=[A;AB];
% b=[b;beq];

AB1=cat(1,AB1,AB);
beq1=cat(1,beq1,beq);


%============Second booster doses for ages 05 to 11 at the end of simulation============%
Cc1=[];
for j=1:N
    Cc=zeros(1,M*O*N);
    for k=1:O
    for i=1:M 
            eta=0.5; 
        Cc((k-1)*N*M+(j-1)*M+i)=1/(2*eta);
    end
    end
    Cc1=cat(1,Cc1,Cc);
end


ABB=[zeros(17,MM(1)*O*N*3) zeros(17,MM(2)*O*N*3) 0*Cc1 0*Cc1 Cc1];
beqbb=Pop_No_VaxxBB_05;
%====================================================================================%
ABB1=cat(1,ABB1,ABB);
beqq1=cat(1,beqq1,beqbb);


AB1=cat(1,AB1,ABB1);
beq1=cat(1,beq1,beqq1);



%% Constraint on Eligible for Booster
%The number of fully vaccinated individuals 12 weeks prior of the vaccination week
%is set as the upper limit.
%===================================For 18 onwards===================================%
C3=[];
b3=[];
M=5;
[~, ~, raw]=xlsread('Constraint_Data_Booster_Jan20.xlsx');
[row, ~]=size(raw);
Raw=raw(17:28,4:20); Raw=cell2mat(Raw); Raw=round(Raw);
VaxxB_Reg_18=raw(17:33,2); VaxxB_Reg_18=cell2mat(VaxxB_Reg_18); 

%For the first 12 weeks of the simulation
for k=1:12
Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;

C3=cat(1, C3, [0*Cbb Cbb 0*Cbb zeros(N,MM(2)*O*N*3) zeros(N,MM(3)*O*N*3)]);
b3=cat(1, b3, Raw(k,:)'-VaxxB_Reg_18); 
%b3=[b3;Raw(k,:)'-VaxxB_Reg_18];
end
%------------------------------------------------

%For the 13th until the last week of the simulation
B3_18=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(1));   
     CC1_18=CC1_18(1:N,1:N*MM(1));
     CC1_18=-CC1_18;
     
     
     for i=1:O-25
        CC1_18v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_18v2=cat(2,CC1_18v2,CC1_18);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end

         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_18v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_18v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
       
        
  C3=cat(1,C3, [CbbP CbbB 0*CbbB zeros(N,MM(2)*O*N*3) zeros(N,MM(3)*O*N*3)]);
  b3=cat(1,b3, B3_18);
        
     end
%----------------------------------------------------------

%====================================================================================%
A=[A;C3];
b=[b;b3];



%====================================For 12 to 17====================================%
C3=[];
b3=[];
M=3;
Raw=raw(34:45,4:20); Raw=cell2mat(Raw); Raw=round(Raw);
VaxxB_Reg_12=raw(17:33,3); VaxxB_Reg_12=cell2mat(VaxxB_Reg_12); 


%For the first 12 weeks of the simulation
for k=1:12
 Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;


C3=cat(1, C3, [zeros(N, MM(1)*O*N*3) 0*Cbb Cbb 0*Cbb zeros(N,MM(3)*O*N*3)]);
b3=cat(1, b3, Raw(k,:)'-VaxxB_Reg_12); 

end
%----------------------------------------------------------

%For the 13th until the last week of the simulation
B3_12=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(2));   
     CC1_12=CC1_12(1:N,1:N*MM(2));
     CC1_12=-CC1_12;
     
     
     for i=1:O-25
        CC1_12v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_12v2=cat(2,CC1_12v2,CC1_12);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end


         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_12v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_12v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
       
        
  C3=cat(1,C3, [zeros(N, MM(1)*O*N*3) CbbP CbbB 0*CbbB zeros(N,MM(3)*O*N*3)]);
  b3=cat(1,b3, B3_12);
         
     end
%-----------------------------------------------------------------------



%====================================================================================%
A=[A;C3];
b=[b;b3];


%====================================For 05 to 11====================================%
C3=[];
b3=[];
M=3;
Raw=raw(48:59,4:20); Raw=cell2mat(Raw); Raw=round(Raw);

%For the first 12 weeks of the simulation
for k=1:12
 Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;

C3=cat(1, C3, [zeros(N, MM(1)*O*N*3) zeros(N,MM(2)*O*N*3) 0*Cbb Cbb 0*Cbb]);
b3=cat(1, b3, Raw(k,:)'-zeros(17,1)); 
end
%----------------------------------------------------------

%For the 13th until the last week of the simulation
B3_05=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(3));   
     CC1_05=CC1_05(1:N,1:N*MM(3));
     CC1_05=-CC1_05;
     
     
     for i=1:O-25
        CC1_05v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_05v2=cat(2,CC1_05v2,CC1_05);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end


         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_05v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_05v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
       
        
  C3=cat(1,C3, [zeros(N, MM(1)*O*N*3) zeros(N,MM(2)*O*N*3) CbbP CbbB 0*CbbB]);
  b3=cat(1,b3, B3_05);

     end
%-----------------------------------------------------------------------
%====================================================================================%

A=[A;C3];
b=[b;b3];



%% Constraint on Eligible for Second Booster
%The number of individuals that have taken their first booster 12 weeks prior of the vaccination week
%is set as the upper limit.
%===================================For 18 onwards===================================%
C3=[];
b3=[];
M=5;
[~, ~, raw]=xlsread('Constraint_Data_SecondBooster_Jan20.xlsx');
[row, ~]=size(raw);
Raw=raw(1:12,4:20); Raw=cell2mat(Raw); Raw=round(Raw);
VaxxB_Reg_18=raw(1:17,1); VaxxB_Reg_18=cell2mat(VaxxB_Reg_18); 

%For the first 12 weeks of the simulation
for k=1:12
Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;

C3=cat(1, C3, [0*Cbb 0*Cbb Cbb zeros(N,MM(2)*O*N*3) zeros(N,MM(3)*O*N*3)]);
b3=cat(1, b3, Raw(k,:)'-VaxxB_Reg_18); 
%b3=[b3;Raw(k,:)'-VaxxB_Reg_18];
end
%------------------------------------------------

%For the 13th until the last week of the simulation
B3_18=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(1));   
     CC1_18=CC1_18(1:N,1:N*MM(1));
     CC1_18=-CC1_18;
     
     
     for i=1:O-13
        CC1_18v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_18v2=cat(2,CC1_18v2,CC1_18);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end

         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_18v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_18v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
       CbbB1=logical(CbbP); CbbB1= double(CbbB1); CbbB1=-CbbB1;
       CbbBB=CbbB; 
       
  C3=cat(1,C3, [0*CbbP CbbB1 CbbBB zeros(N,MM(2)*O*N*3) zeros(N,MM(3)*O*N*3)]);
  b3=cat(1,b3, B3_18);
        
     end
%----------------------------------------------------------

%====================================================================================%
A=[A;C3];
b=[b;b3];

%[activeA,A,b,AB1,beq1] = iisremover(coeff,A,b,AB1,beq1,zeros(length(coeff),1),max(Vsup)+zeros(length(coeff),1));
%====================================For 12 to 17====================================%
C3=[];
b3=[];
M=3;
Raw=raw(18:29,4:20); Raw=cell2mat(Raw); Raw=round(Raw);
VaxxB_Reg_12=raw(18:34,2); VaxxB_Reg_12=cell2mat(VaxxB_Reg_12); 


%For the first 12 weeks of the simulation
for k=1:12
 Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;


C3=cat(1, C3, [zeros(N, MM(1)*O*N*3) 0*Cbb 0*Cbb Cbb zeros(N,MM(3)*O*N*3)]);
b3=cat(1, b3, Raw(k,:)'-VaxxB_Reg_12); 

end
%----------------------------------------------------------

%For the 13th until the last week of the simulation
B3_12=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(2));   
     CC1_12=CC1_12(1:N,1:N*MM(2));
     CC1_12=-CC1_12;
     
     
     for i=1:O-13
        CC1_12v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_12v2=cat(2,CC1_12v2,CC1_12);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end


         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_12v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_12v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
       CbbB1=logical(CbbP); CbbB1= double(CbbB1); CbbB1=-CbbB1;
       CbbBB=CbbB;
        
  C3=cat(1,C3, [zeros(N, MM(1)*O*N*3) 0*CbbP CbbB1 CbbBB zeros(N,MM(3)*O*N*3)]);
  b3=cat(1,b3, B3_12);
         
     end
%-----------------------------------------------------------------------



%====================================================================================%
A=[A;C3];
b=[b;b3];


%====================================For 05 to 11====================================%
C3=[];
b3=[];
M=3;
Raw=raw(32:43,4:20); Raw=cell2mat(Raw); Raw=round(Raw);

%For the first 12 weeks of the simulation
for k=1:12
 Cbb=[];
  
   for j=1:N
     if k==1 
            Cb=zeros(1,M*O*N);
     else
            Cb=Cbb2(j,:);
     end  
    for i=1:M
        Cb((k-1)*N*M+(j-1)*M+i)=1;
    end   
 Cbb=cat(1,Cbb,Cb);
   end
Cbb2=Cbb;

C3=cat(1, C3, [zeros(N, MM(1)*O*N*3) zeros(N,MM(2)*O*N*3) 0*Cbb 0*Cbb Cbb]);
b3=cat(1, b3, Raw(k,:)'-zeros(17,1)); 
end
%----------------------------------------------------------

%For the 13th until the last week of the simulation
B3_05=b3(188:204);
     CBB2=Cbb2(1:N,1:N*MM(3));   
     CC1_05=CC1_05(1:N,1:N*MM(3));
     CC1_05=-CC1_05;
     
     
     for i=1:O-13
        CC1_05v2=[]; 
        CBB2v2=[];
         for ii=1:i
            CC1_05v2=cat(2,CC1_05v2,CC1_05);
            CBB2v2=cat(2,CBB2v2,CBB2);
         end


         CBB2v3=1*CBB2v2;
         [~,colB]=size(CBB2v3);
         [~, colP]=size(CC1_05v2);
        CbbP=0*Cbb; CbbP(:,1:colP)=[];
        CbbP=cat(2,CC1_05v2,CbbP); 
      
        CbbB=Cbb2;
        
        CbbB(:,N*M*12+1:N*M*12+colB)=CBB2v2;   
     
        CbbB(:,1:colB)=[];
        CbbB=cat(2,CBB2v3,CbbB);
        
        CbbB1=logical(CbbP); CbbB1= double(CbbB1); CbbB1=-CbbB1;
        CbbBB=CbbB;     
        
  C3=cat(1,C3, [zeros(N, MM(1)*O*N*3) zeros(N,MM(2)*O*N*3) 0*CbbP CbbB1 CbbBB]);
  b3=cat(1,b3, B3_05);

     end
%-----------------------------------------------------------------------
%====================================================================================%

A=[A;C3];
b=[b;b3];





%% Constraint on Max Weekly Vaccination
%Set to at most 5,040,000 of inoculated vaccine doses each week of the simulation.
%===================================For 18 onwards===================================%
M=5;
AA=[];
bb=[];
for k=1:O
    C4=zeros(1,M*O*N);
    C4b=C4;
    C4bb=C4; 
    
    for j=1:N
        for i=1:M
            if k<=O-24
        C4((k-1)*N*M+(j-1)*M+i)=1;
            end
            
            if k<=O-12
        C4b((k-1)*N*M+(j-1)*M+i)=1;
            end
        C4bb((k-1)*N*M+(j-1)*M+i)=1;
        end
    end
    AA=cat(1,AA, [C4 C4b C4bb]);
    bb=cat(1,bb, 0); 
end
AA1=AA;
%====================================================================================%


%====================================For 12 to 17====================================%
M=3;
AA=[];

for k=1:O
    C4=zeros(1,M*O*N);
    C4b=C4;
    C4bb=C4;
    
    for j=1:N
        for i=1:M
            if k<=O-24
        C4((k-1)*N*M+(j-1)*M+i)=1;
            end
            
            if k<=O-12
        C4b((k-1)*N*M+(j-1)*M+i)=1;
            end
        C4bb((k-1)*N*M+(j-1)*M+i)=1;
        
        end
    end
    AA=cat(1,AA, [C4 C4b C4bb]);

end
AA1=cat(2,AA1,AA);
%====================================================================================%


%====================================For 05 to 11====================================%
M=3;
eta=1;
AA=[];

for k=1:O
    C4=zeros(1,M*O*N);
    C4b=C4;
    C4bb=C4;
         
    for j=1:N
        for i=1:M
               if k<=O-24
        C4((k-1)*N*M+(j-1)*M+i)=1;
               end
               
               if k<=O-12
        C4b((k-1)*N*M+(j-1)*M+i)=1;
               end
        C4bb((k-1)*N*M+(j-1)*M+i)=1;
        end
    end
    AA=cat(1,AA, [C4 C4b C4bb]);

end
AA1=cat(2,AA1,AA);
%====================================================================================%
bb1=bb*0+5040000; %Weekly Vaccine Dose Limit 
A=cat(1,A,AA1);
b=cat(1,b,bb1);


%% Constraint on Same-Second Dose for Primary Vaccination
%Pfizer, Moderna, Astrazeneca, and Sinovac vaccines are the covid vaccine brands
%that needed two shots (or doses) to get fully vaccinated. Pfizer vaccine
%is set to receive the second dose 3 weeks after the first dose. The other
%three brands are set to take the second dose 4 weeks apart from the first dose.
%===================================For 18 onwards===================================%
A_samedose=[];
b_samedose=[];
M=5;
for k=1:O-27
    i=1; 
       for j=1:N
    C4=zeros(1,M*O*N);
    C5=C4;
    C4((k-1)*N*M+(j-1)*M+i)=1;
    C4((k+2)*N*M+(j-1)*M+i)=-1;
    A_samedose=cat(1,A_samedose, [C4 0*C4 0*C4 zeros(1,MM(2)*O*N*3) zeros(1,MM(3)*O*N*3)]);
    b_samedose=cat(1,b_samedose, 0);
       end
end

A=[A;A_samedose];
b=[b;b_samedose];

vax4=[2;4;5];

A_samedose2=[];
b_samedose2=[];


for k=1:O-28

    for iii=1:3
    i=vax4(iii);
        for j=1:N
        C4=zeros(1,M*O*N);
        C5=C4;
        C4((k-1)*N*M+(j-1)*M+i)=1;
        C4((k+3)*N*M+(j-1)*M+i)=-1;
        A_samedose2=cat(1,A_samedose2, [C4 0*C4 0*C4 zeros(1,MM(2)*O*N*3) zeros(1,MM(3)*O*N*3)]);
        b_samedose2=cat(1,b_samedose2,0);

        end
    end
end
A=[A;A_samedose2];
b=[b;b_samedose2];
%====================================================================================%

%====================================For 12 to 17====================================%
A_samedose=[];
b_samedose=[];
M=3;
for k=1:O-27
    i=1; 
       for j=1:N
    C4=zeros(1,M*O*N);
    C5=C4;
    C4((k-1)*N*M+(j-1)*M+i)=1;
    C4((k+2)*N*M+(j-1)*M+i)=-1;
    A_samedose=cat(1,A_samedose, [zeros(1,MM(1)*O*N*3) C4 0*C4 0*C4 zeros(1,MM(3)*O*N*3)]);
    b_samedose=cat(1,b_samedose, 0);
       end
end

A=[A;A_samedose];
b=[b;b_samedose];

vax4=[2;3];

A_samedose2=[];
b_samedose2=[];


for k=1:O-28

    for iii=1:2
        
    i=vax4(iii);
        for j=1:N
        C4=zeros(1,M*O*N);
        C5=C4;
        C4((k-1)*N*M+(j-1)*M+i)=1;
        C4((k+3)*N*M+(j-1)*M+i)=-1;
        A_samedose2=cat(1,A_samedose2, [zeros(1,MM(1)*O*N*3) C4 0*C4 0*C4 zeros(1,MM(3)*O*N*3)]);
        b_samedose2=cat(1,b_samedose2,0);

        end
    end
end
A=[A;A_samedose2];
b=[b;b_samedose2];
%====================================================================================%

%====================================For 05 to 11====================================%
A_samedose=[];
b_samedose=[];
M=3;
for k=1:O-27

    i=1; 
       for j=1:N
    C4=zeros(1,M*O*N);
    C5=C4;
    C4((k-1)*N*M+(j-1)*M+i)=1;
    C4((k+2)*N*M+(j-1)*M+i)=-1;
    A_samedose=cat(1,A_samedose, [zeros(1,MM(1)*O*N*3) zeros(1,MM(2)*O*N*3) C4 0*C4 0*C4]);
    b_samedose=cat(1,b_samedose, 0);

       end
end

A=[A;A_samedose];
b=[b;b_samedose];

vax4=2;

A_samedose2=[];
b_samedose2=[];


for k=1:O-28

    for iii=1:1
    i=vax4(iii);
        for j=1:N
        C4=zeros(1,M*O*N);
        C5=C4;
        C4((k-1)*N*M+(j-1)*M+i)=1;
        C4((k+3)*N*M+(j-1)*M+i)=-1;
        A_samedose2=cat(1,A_samedose2, [zeros(1,MM(1)*O*N*3) zeros(1,MM(2)*O*N*3) C4 0*C4 0*C4]);
        b_samedose2=cat(1,b_samedose2,0);

        end
    end
end
A=[A;A_samedose2];
b=[b;b_samedose2];
%====================================================================================%




%%
%SOLVING

options = optimoptions('linprog','Algorithm','interior-point','MaxIterations',100*(size(AB1,1)+size(A,1)+length(coeff)),'OptimalityTolerance',1e-10); %interior-point algorithm is used to obtain the solution
X=linprog(coeff,A,b,AB1,beq1,zeros(length(coeff),1),max(Vsup)+zeros(length(coeff),1),options); X=round(X);
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------
toc;
tic;
%%
%Display Results
Region_list={'CAR','I','II','III','IVA','NCR','IVB','V','VI','VII','VIII','IX','X','XI','XII','XIII','BARMM'};
%===================================For 18 onwards===================================%
M=5;
Y=reshape(X(1:MM(1)*O*N),MM(1)*N,O);
Z=reshape(X(MM(1)*O*N+1:2*MM(1)*O*N),MM(1)*N,O);
ZZ=reshape(X(2*MM(1)*O*N+1:3*MM(1)*O*N),MM(1)*N,O);
titlelist={'Pfizer','Moderna','Janssen', 'Astrazeneca', 'Sinovac'};
fig=figure;
for i=1:M*N
    subplot(17,M,i);
    plot(Y(i,:));
    hold on;
    subplot(17,M,i)
    plot(Z(i,:));
    subplot(17,M,i)
    plot(ZZ(i,:));
    hold on;
    if i>=1 && i<6
    title(titlelist(i));
    end
    xlim([0 53]);
    if mod(i,M)==1
       ylabel(Region_list(ceil(i/M))) 
    end
end
suptitle('For Ages 18 and above');
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
Ydaw=ylabel(han,'REGIONS','fontweight','bold','fontsize', 12,'color','k');
set(Ydaw, 'Units', 'Normalized', 'Position', [-0.025, 0.5, 1]);
%====================================================================================%
Y1=Y;
Z1=Z;
ZZ1=ZZ;

%====================================For 12 to 17====================================%
M=3;
Y=reshape(X(MM(1)*O*N*3+1:MM(1)*O*N*3+MM(2)*O*N),M*N,O);
Z=reshape(X(MM(1)*O*N*3+MM(2)*O*N+1:MM(1)*O*N*3+MM(2)*O*N*2),M*N,O);
ZZ=reshape(X(MM(1)*O*N*3+MM(2)*O*N*2+1:MM(1)*O*N*3+MM(2)*O*N*3),M*N,O);
titlelist={'Pfizer','Moderna', 'Sinovac'};
fig1=figure;
for i=1:M*N
    subplot(17,M,i);
    plot(Y(i,:));
    hold on;
    subplot(17,M,i)
    plot(Z(i,:));
    hold on
    subplot(17,M,i)
    plot(ZZ(i,:));    
    hold on;
     if i>=1 && i<4
    title(titlelist(i));
    end
    xlim([0 53]);
    if mod(i,M)==1
       ylabel(Region_list(ceil(i/M))) 
    end   
end
suptitle('For Ages 12 to 17');
han=axes(fig1,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
Ydaw=ylabel(han,'REGIONS','fontweight','bold','fontsize', 12,'color','k');
set(Ydaw, 'Units', 'Normalized', 'Position', [-0.025, 0.5, 1]);
%====================================================================================%
Y2=Y;
Z2=Z;
ZZ2=ZZ;

%====================================For 05 to 11====================================%
M=3;
Y=reshape(X(MM(1)*O*N*3+MM(2)*O*N*3+1:MM(1)*O*N*3+MM(2)*O*N*3+MM(3)*O*N),M*N,O);
Z=reshape(X(MM(1)*O*N*3+MM(2)*O*N*3+MM(3)*O*N+1:MM(1)*O*N*3+MM(2)*O*N*3+MM(3)*O*N*2),M*N,O);
ZZ=reshape(X(MM(1)*O*N*3+MM(2)*O*N*3+MM(3)*O*N*2+1:length(X)),M*N,O);
titlelist={'Pfizer','Moderna', 'Sinovac'};
fig2=figure;
for i=1:M*N
    subplot(17,M,i);
    plot(Y(i,:));
    hold on;
    subplot(17,M,i)
    plot(Z(i,:));
    hold on;
    subplot(17,M,i)
    plot(ZZ(i,:));    
    hold on;    
     if i>=1 && i<4
    title(titlelist(i));
    end
    xlim([0 53]);
    if mod(i,M)==1
       ylabel(Region_list(ceil(i/M))) 
    end    
end
suptitle('For Ages 05 to 11');
han=axes(fig2,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
Ydaw=ylabel(han,'REGIONS','fontweight','bold','fontsize', 12,'color','k');
set(Ydaw, 'Units', 'Normalized', 'Position', [-0.025, 0.5, 1]);
%====================================================================================%
Y3=Y;
Z3=Z;
ZZ3=ZZ;



toc;

[results,sum_tanan]=copy_save_results(Y1,Y2,Y3,Z1,Z2,Z3,ZZ1,ZZ2,ZZ3);

disp(sum_tanan); % displays (in order the) [overall doses needed; total pfizer doses administered;
                 %                         total moderna doses given; total Janssen doses used; 
                 %                         total astrazeneca/sputnik V shots; total sinovac doses inoculated]
                 % to minimize the infection cause by the Covid-19 disease through vaccinating 
                 % all the eligible population in Philippine Regions  

