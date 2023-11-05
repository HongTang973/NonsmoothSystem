% construct a 3D example to show the parameter space with specific
% features

% ------SubCASE 2: LAMBDA_3 = -1;  lambda_1 & lambda_2 are complex -------%

%*************** start parallel computation ***************%
% find the numbers of this
Initialization_for_Par;


% give span of three parameteres and mesh grid
d_mesh = 0.1;
X_span = [-1:d_mesh:1];
Y_span = [d_mesh:d_mesh:1]; % w part
Z_span = [-1:d_mesh:1];

% discard zero values
X_span(X_span==0) = [];
Y_span(Y_span==0) = [];
Z_span(Z_span==0) = [];
% two ways to produce the sampling points: random / uniform grid

[X,Y,Z] = meshgrid(X_span,Y_span,Z_span);

% for every point of parameter: try with the algorithm and record the
% results

para_list = [reshape(X,[],1), reshape(Y,[],1), reshape(Z,[],1)];



lam1 = para_list(:,1)+ i*para_list(:,2);
lam2 = para_list(:,1)- i*para_list(:,2);
b_3  = para_list(:,3);
% preallocate a table to store the results.


varNames = {'lambda_1','lambda_2','lambda_3','b2','b3','num of LCO', ...
    'roots','F Value','sign V', ...
    'Max', 'datetime'};

varTypes = {'double','double','double','double','double','double',...
    'double', 'double','double', 'double', 'datetime'};
sz = [length(para_list) length(varNames)];
Monte_3D = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);


tic
parfor i = 1:length(para_list)
    equi_type = 1;
    lambda_1 = lam1(i);
    lambda_2 = lam2(i);
    lambda_3 = -1;
    b2 = 1.75;
    b3 = b_3(i);
% % for i =1
%     equi_type = 1;
%     lambda_1 = 0.9;
%     lambda_2 = 0.6;
%     lambda_3 = -1;
%     b2 = 1.75;
%     b3 = 0.2;
    C = [1,0,0];
    
%     try
        [A,R,C]= Matrices_3D_impact(lambda_1,lambda_2,lambda_3,b2,b3,C);
        
        % set the searching limit
        [V1,D1]=eig(A);
        d = real(diag(D1));
        w = abs(imag(diag(D1)));
        
        % due to the all real eigenvalues, no w
        UT = pi;
        % determine the sampling frequency
        Omega = max(abs(diag(D1)));
        fs = max( 10*ceil(2*(Omega/2/pi)), 1000);
        
        % for different Evolution time
        a =0;
        b = UT;
        delta = 1/fs;
        T = 1/fs:delta:b;
        %
        MAX = zeros(1,length(T));
        F_1 = zeros(1,length(T));
        V_sign = zeros(1,length(T));
        LOCI= zeros(size(A,1),length(T));
        vector= zeros(size(A,1),length(T));
        
        % began searching
        for j=1:length(T)
            [V_sign(j),LOCI(:,j),MAX(j),vector(:,j),F_1(j),~] = LCO_Det_search(T(j),R,A,C,equi_type);
        end
        index0 =sign(F_1);
        index1 = abs(diff(index0))>0;
        % filter the singularity case
        index_s = abs(diff(F_1))/delta < (1/delta);
        index1 =index1 & index_s;
        index2 = [index1,0];
        index3 = [0,index1];
        index2=find(index2==1);
        index3=find(index3==1);
        % section partia
        ratio = abs(F_1(index2))./(abs(F_1(index2))+abs(F_1(index3)));
        T_chosen =(1-ratio).*T(index2)+ratio.*T(index3);
        MAX_chosen = (1-ratio).*MAX(index2)+ratio.*MAX(index3);
        F1_chosen = (1-ratio).*F_1(index2)+ratio.*F_1(index3);
        Sign_chosen =(1-ratio).*V_sign(index2)+ratio.*V_sign(index3);
        
        index4 =  (abs(F1_chosen)<1e-3) & (T_chosen>1e-3);
        
        T_chosen  =T_chosen(index4);
        F1_chosen =F1_chosen(index4);
        
        if ~isempty(T_chosen)
            % Number of roots found:
            % the roots
            % the sign of the velocity related to the LCO (to distinguish to SAE or SPE type)
            % the stability of the root
            V_sign_LCO = [];
            MAX = [];
            LCO = [];
            F_value = [];
            
            for k=1:length(T_chosen)
                % get the LCO with statevariables
                [op1,~,~,op4,op5,~] = LCO_Det_search(T_chosen(k),R,A,C,equi_type);
                
                %check the floque multipliers of the foud LCO
                [Mono_p,Salt_p]=Floque_Multipliers(T_chosen(k),op4,R,A,C);
                
                V_sign_LCO = [V_sign_LCO,op1];
                MAX = [MAX,max(abs(Salt_p)) ];
                LCO  = [ LCO; op4];
                F_value = [F_value,op5];
            end
            
            Monte_3D(i,:) = {lambda_1,lambda_2,lambda_3,b2,b3, length(T_chosen), ...
                T_chosen,F_value, V_sign_LCO, MAX , datetime};
            
        else
            
            Monte_3D(i,:) = { lambda_1,lambda_2,lambda_3,b2,b3, 0, ...
                0 ,1, 2, -5 , datetime};
            
        end
        %
%     catch ME
%         
%         disp(['error at',num2str(lambda_1),',',num2str(lambda_2),',',num2str(b3)])
%         
%     end
end

CPU_tiem = toc

save('Monte_3D_case2.mat','Monte_3D')

