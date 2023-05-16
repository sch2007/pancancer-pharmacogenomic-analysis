%% Initialize
% replace Init with addpath to the LOBICO directory
Init

disp('***** Running LOBICO *****');
extraLabel = '';

drugList = {'Afatinib','Buparlisib','Cobimetinib','Crizotinib',...
            'Dabrafenib', 'Dinaciclib', 'Erlotinib', 'Everolimus',...
            'Apitolisib', 'Ganetespib', 'Gedatolisib', 'LY2603618',...
            'Lapatinib', 'Sapanisertib', 'Paclitaxel', 'Temsirolimus',...
            'Trametinib', 'Vemurafenib', 'Vorinostat'};

drugFilename = 'Supplementary Table 5 LOBICO data for paper.xlsx';

extraLabel = sprintf('%s_t0.1', extraLabel);

verbosity = 1;  % 0..no output, 1..minimal output, 2..full output

% -- default --
log_units = true;
normalized = false;

if normalized
    extraLabel = sprintf('%s_norm', extraLabel);
end
if log_units
    extraLabel = sprintf('%s_log10', extraLabel);
end
thresholdFilename = sprintf('AUC%s.csv', extraLabel);
disp(sprintf('* AUC: log10 = %d, normalized = %d', log_units, normalized))

%% Logic model complexity

% first entry is K, second is M
KMs = [[1, 1]; [1, 2]; [2, 1]; [1, 3]; [3, 1]; [2, 2]; [1, 4]; [4, 1]];
extraLabel = sprintf('%s_limit4', extraLabel);

KM_max = max(max(KMs));

% sens & spec run from 0 to 1 in step sizes of d_sens, d_spec
d_sens = 0.05;
d_spec = 0.05;

% No of folds in cross validation for hyperparameter tuning (K, M)
% i.e., split of training data into training(proper) and validation
nfold = 4;

% No of folds in cross validations for testing (training-test split)
niter = 5;

% load threshold value
fileID = fopen(thresholdFilename);
thresholdData = textscan(fileID,'%s%d%f%*[^\n]', 'CommentStyle','#', 'Delimiter', ',');
fclose(fileID);

outFilename1 = sprintf('Validation_n%d_k%d_KM%d%s.csv', niter, nfold, KM_max, extraLabel);
outfile1 = fopen(outFilename1, 'a');
outFilename2 = sprintf('Testing_n%d_KM%d%s.csv', niter, KM_max, extraLabel);
outfile2 = fopen(outFilename2, 'a');

fprintf(outfile1, "drug,cell_lines,genes,lg_thresh,niter,iiter,clauses,literals,min_sens,min_spec,nfold,ifold," + ...
                  "num_true_pos_train,num_pos_train,num_true_neg_train,num_neg_train,Matthews_train,Kappa_train,F1_train," + ...
                  "num_true_pos_val,num_pos_val,num_true_neg_val,num_neg_val,Matthews_val,Kappa_val,F1_val,formula\n");

fprintf(outfile2, "drug,cell_lines,genes,lg_thresh,niter,iiter,clauses,literals,min_sens,min_spec," + ...
                  "num_true_pos_trainval,num_pos_trainval,num_true_neg_trainval,num_neg_trainval,Matthews_trainval,Kappa_trainval,F1_trainval," + ...
                  "num_true_pos_test,num_pos_test,num_true_neg_test,num_neg_test,Matthews_test,Kappa_test,F1_test,formula\n");

%% Outer loop over drugs (completely independent for each drug)
for iDrug = 1:length(drugList)
    
    drug = drugList{iDrug};
    
    if verbosity > 0
        disp(sprintf('**** Drug %s ****', drug));
    end   

    for iidrug = 1:length(thresholdData{1})
        drug2 = strrep(thresholdData{1}{iidrug}, '"', '');
        if strcmp(drug, drug2)
           th = thresholdData{3}(iidrug);
        end
    end
    if verbosity > 0
        disp(sprintf('* Binarization threshold %g *', th));
    end

    %% Load data
    [Samples, Features, AUCs, MutationMatrix] = ...
        parsexls2(drugFilename, drug, 1, 3);
    
    if normalized
        AUCs = AUCs / max(AUCs);
    end
    if log_units
        AUCs = log10(AUCs);
    end

    % remove cell-lines that do not have complete set of gene features
    missing_list = [];
    for iSample = 1:length(Samples)
        if any(ismissing(MutationMatrix(iSample, :)) | ismissing(AUCs(iSample)))
            missing_list = [missing_list, iSample];
        end
    end
    MutationMatrix(missing_list, :) = [];
    Samples(missing_list) = [];
    AUCs(missing_list) = [];
    if verbosity > 0
        disp(sprintf('* # cell lines = %d *', length(Samples)))
        disp(sprintf('* # gene features = %d *', length(Features)))
    end

    %% Cross validation
    rng(1); % so that cross-validation can be reproduced
    [testing_indicess, trainval_indicess] = balanced_crossval(length(Samples), niter, AUCs < th);
    for iiter = 1:niter  % cross-validation for training-test split
        if verbosity > 0
            disp(sprintf('* CV iteration: %d / %d', iiter, niter))
        end
        
        testing_indices = testing_indicess{iiter};
        trainval_indices = trainval_indicess{iiter};

        [validation_indicess, training_indicess] = balanced_crossval(length(trainval_indices), nfold, AUCs(trainval_indices) < th);

        for iKM = 1:length(KMs)  % model complexity

            K = KMs(iKM, 1);
            M = KMs(iKM, 2);
            
            if verbosity > 1
                disp(sprintf('K (#clauses) = %d, M (#literals) = %d', K, M));
            end

            min_specs = 0:d_spec:1;
            for ispec = 1:length(min_specs)
                min_spec = min_specs(ispec);

                min_senss = (1-min_spec):d_sens:1;
                for isens = 1:length(min_senss)
                    min_sens = min_senss(isens);
                    if verbosity > 1
                        disp(sprintf('* Minimum sensitivity = %g *', min_sens));
                        disp(sprintf('* Minimum specificity = %g *', min_spec));
                    end

                    for ifold = 1:nfold  % splitting training into training/validation
                        validation_indices = trainval_indices(validation_indicess{ifold});
                        training_indices = trainval_indices(training_indicess{ifold});

                        %% Create binary input, output and weight vector

                        %binary input
                        X = MutationMatrix(training_indices, :);
                        [N,P] = size(X);

                        %binarization threshold th
                        Y = double(AUCs(training_indices)<th);
                        W = abs(AUCs(training_indices)-th);

                        %class weights
                        FPW = 1;                  %Total weight on positive class (Y==1)
                        FPN = 1;                  %Total weight on negative class (Y==0)

                        %normalize weights
                        W(Y==1) = FPW*W(Y==1)./sum(W(Y==1));
                        W(Y~=1) = -(FPN*W(Y~=1)./sum(W(Y~=1)));

                        %% Cplex options
                        param = cat(1,{'timelimit.Cur',60,'MaxTime'},...                            %Maximum time for IP (in seconds)
                                      {'mip.tolerances.integrality.Cur',1e-5,'Integrality'},...     %Integrality contraint; default = 1e-5 (see cplex.Param.mip.tolerances.integrality.Help)
                                      {'mip.tolerances.mipgap.Cur',1e-4,'RelGap'},...               %Optimality gap tolerance; default = 1e-4 (0.01% of optimal solution, set to 0.05, 5% for early termination, approximate solution) (see cplex.Param.mip.tolerances.mipgap.Help)
                                      {'threads.Cur',8,'Threads'},...                               %Number of threads to use (default = 0, automatic) (see  cplex.Param.threads.Help);
                                      {'parallel.Cur',-1,'ParallelMode'},...                        %Parallel optimization mode,  -1 = opportunistic 0 = automatic 1 = deterministic (see cplex.Param.parallel.Help)
                                      {'mip.pool.relgap.Cur',1e-1,'Pool_relgap'},...                %Relative gap for suboptimal solutions in the pool; default 0.1 (10%)
                                      {'mip.pool.intensity.Cur',1,'Pool_intensity'},...             %Pool intensity; default 1 = mild: generate few solutions quickly (see  cplex.Param.mip.pool.intensity.Help)
                                      {'mip.pool.replace.Cur',2,'Pool_replace'},...                 %Pool replacement strategy; default 2 = replace least diverse solutions (see  cplex.Param.mip.pool.replace.Help)
                                      {'mip.pool.capacity.Cur',11,'Pool_capacity'},...              %Pool capacity; default 11 = best + 10 (see  cplex.Param.mip.pool.replace.Help)
                                      {'mip.limits.populate.Cur',11,'Pool_size'});                  %Number of solutions generated; default 11 = best + 10 (see  cplex.Param.mip.limits.populate.Help)

                        %% Cplex solver
                        sol = lobico(X,W,K,M,1, param, min_spec, min_sens);

                        %% check whether solution found
                        if ~isfield(sol.Solution, 'x')
                            continue
                        end

                        %% Check solution
                        %inferred formula
                        x = round(sol.Solution.x);
                        SolMat = getsolution(x,K,M,P);
                        str = showformula(SolMat,K,M,Features);
                        if isempty(str)
                            str = '0';
                        else
                            str = strtrim(str);
                        end
                        if verbosity > 1
                            disp('***********************');
                            disp('Inferred logic model');
                            disp(str);
                            disp('***********************');
                        end

                        %% Apply model to training data
                        X_train = MutationMatrix(training_indices, :);
                        Y_train = double(AUCs(training_indices)<th);

                        %apply the trained model to the data
                        Y_train_pred = applymodel(x, X_train, K, M, P);

                        num_true_pos_train = sum(Y_train_pred(find(Y_train==1)));
                        num_pos_train = sum(Y_train);

                        num_true_neg_train = sum(1 - Y_train_pred(find(Y_train==0)));
                        num_neg_train = sum(1 - Y_train);

                        [c_matrixp,Result_train]= confusion.getMatrix(~Y_train, ~Y_train_pred, 0);

                        %% Apply model to validation data
                        X_val = MutationMatrix(validation_indices, :);
                        Y_val = double(AUCs(validation_indices)<th);

                        %apply the trained model to the data
                        Y_val_pred = applymodel(x, X_val, K, M, P);

                        num_true_pos_val = sum(Y_val_pred(find(Y_val==1)));
                        num_pos_val = sum(Y_val);

                        num_true_neg_val = sum(1 - Y_val_pred(find(Y_val==0)));
                        num_neg_val = sum(1 - Y_val);

                        [c_matrixp, Result_val]= confusion.getMatrix(~Y_val, ~Y_val_pred, 0);

                        fprintf(outfile1, '%s,%d,%d,%g,%d,%d,%d,%d,%g,%g,%d,%d,%d,%d,%d,%d,%g,%g,%g,%d,%d,%d,%d,%g,%g,%g,"%s"\n',...
                                drug, length(Samples), length(Features), th, niter, iiter, K, M, min_sens, min_spec, nfold, ifold, ...
                                num_true_pos_train, num_pos_train, num_true_neg_train, num_neg_train, ...
                                Result_train.MatthewsCorrelationCoefficient, Result_train.Kappa, Result_train.F1_score, ...
                                num_true_pos_val, num_pos_val, num_true_neg_val, num_neg_val, ...
                                Result_val.MatthewsCorrelationCoefficient, Result_val.Kappa, Result_val.F1_score, ...
                                str);

                    end

                    %% train on full training set
                    %binary input
                    X = MutationMatrix(trainval_indices, :);
                    [N,P] = size(X);

                    %binarization threshold th
                    Y = double(AUCs(trainval_indices)<th);
                    W = abs(AUCs(trainval_indices)-th);

                    %class weights
                    FPW = 1;                  %Total weight on positive class (Y==1)
                    FPN = 1;                  %Total weight on negative class (Y==0)

                    %normalize weights
                    W(Y==1) = FPW*W(Y==1)./sum(W(Y==1));
                    W(Y~=1) = -(FPN*W(Y~=1)./sum(W(Y~=1)));

                    %% Cplex options
                    param = cat(1,{'timelimit.Cur',60,'MaxTime'},...                            %Maximum time for IP (in seconds)
                                  {'mip.tolerances.integrality.Cur',1e-5,'Integrality'},...     %Integrality contraint; default = 1e-5 (see cplex.Param.mip.tolerances.integrality.Help)
                                  {'mip.tolerances.mipgap.Cur',1e-4,'RelGap'},...               %Optimality gap tolerance; default = 1e-4 (0.01% of optimal solution, set to 0.05, 5% for early termination, approximate solution) (see cplex.Param.mip.tolerances.mipgap.Help)
                                  {'threads.Cur',8,'Threads'},...                               %Number of threads to use (default = 0, automatic) (see  cplex.Param.threads.Help);
                                  {'parallel.Cur',-1,'ParallelMode'},...                        %Parallel optimization mode,  -1 = opportunistic 0 = automatic 1 = deterministic (see cplex.Param.parallel.Help)
                                  {'mip.pool.relgap.Cur',1e-1,'Pool_relgap'},...                %Relative gap for suboptimal solutions in the pool; default 0.1 (10%)
                                  {'mip.pool.intensity.Cur',1,'Pool_intensity'},...             %Pool intensity; default 1 = mild: generate few solutions quickly (see  cplex.Param.mip.pool.intensity.Help)
                                  {'mip.pool.replace.Cur',2,'Pool_replace'},...                 %Pool replacement strategy; default 2 = replace least diverse solutions (see  cplex.Param.mip.pool.replace.Help)
                                  {'mip.pool.capacity.Cur',11,'Pool_capacity'},...              %Pool capacity; default 11 = best + 10 (see  cplex.Param.mip.pool.replace.Help)
                                  {'mip.limits.populate.Cur',11,'Pool_size'});                  %Number of solutions generated; default 11 = best + 10 (see  cplex.Param.mip.limits.populate.Help)

                    %% Cplex solver
                    sol = lobico(X,W,K,M,1, param, min_spec, min_sens);

                    %% check whether solution found
                    if ~isfield(sol.Solution, 'x')
                        continue
                    end

                    %% Check solution
                    %inferred formula
                    x = round(sol.Solution.x);
                    SolMat = getsolution(x,K,M,P);
                    str = showformula(SolMat,K,M,Features);
                    if isempty(str)
                        str = '0';
                    else
                        str = strtrim(str);
                    end
                    if verbosity > 1
                        disp('***********************');
                        disp('Inferred logic model');
                        disp(str);
                        disp('***********************');
                    end

                    %% Apply model to training+validation data
                    X_trainval = MutationMatrix(trainval_indices, :);
                    Y_trainval = double(AUCs(trainval_indices)<th);

                    %apply the trained model to the data
                    Y_trainval_pred = applymodel(x, X_trainval, K, M, P);

                    num_true_pos_trainval = sum(Y_trainval_pred(find(Y_trainval==1)));
                    num_pos_trainval = sum(Y_trainval);

                    num_true_neg_trainval = sum(1 - Y_trainval_pred(find(Y_trainval==0)));
                    num_neg_trainval = sum(1 - Y_trainval);

                    [c_matrixp,Result_trainval]= confusion.getMatrix(~Y_trainval, ~Y_trainval_pred, 0);

                    %% Apply model to test data
                    X_test = MutationMatrix(testing_indices, :);
                    Y_test = double(AUCs(testing_indices)<th);

                    %apply the trained model to the data
                    Y_pred = applymodel(x, X_test, K, M, P);

                    num_true_pos = sum(Y_pred(find(Y_test==1)));
                    num_pos = sum(Y_test);

                    num_true_neg = sum(1 - Y_pred(find(Y_test==0)));
                    num_neg = sum(1 - Y_test);

                    [c_matrixp,Result]= confusion.getMatrix(~Y_test, ~Y_pred, 0);

                    fprintf(outfile2, '%s,%d,%d,%g,%d,%d,%d,%d,%g,%g,%d,%d,%d,%d,%g,%g,%g,%d,%d,%d,%d,%g,%g,%g,"%s"\n',...
                        drug, length(Samples), length(Features), th, niter, iiter, K, M, min_sens, min_spec, ...
                        num_true_pos_trainval, num_pos_trainval, num_true_neg_trainval, num_neg_trainval, ...
                        Result_trainval.MatthewsCorrelationCoefficient, Result_trainval.Kappa, Result_trainval.F1_score, ...
                        num_true_pos, num_pos, num_true_neg, num_neg, ...
                        Result.MatthewsCorrelationCoefficient, Result.Kappa, Result.F1_score, ...
                        str);

                end
            end
        end
    end
end

fclose(outfile1);
fclose(outfile2);
