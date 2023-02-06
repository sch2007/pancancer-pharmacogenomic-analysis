function [testing_indices, training_indices] = balanced_crossval(n, k, cond)
% balanced_crossval Split x into k folds approx. balance acc. to cond
% i.e., if you have two types of data (cond=0 & cond=1), then
% this code will try to give you folds that spread both types of data
% about equally over the folds

    x = 1:n;

    perm = randperm(n);
    cond = cond(perm);
    
    indices_0 = find(cond==0);
    indices_1 = find(cond==1);

    testing_indices = {};
    training_indices = {};

    for ifold = 1:k

        [testing_0, training_0] = get_indices(ifold, k, sum(cond==0), 0);
        [testing_1, training_1] = get_indices(ifold, k, sum(cond==1), 1);
       
        ind = zeros(n, 1);
        ind(indices_0(testing_0)) = 0;
        ind(indices_0(training_0)) = 1;
        ind(indices_1(testing_1)) = 0;
        ind(indices_1(training_1)) = 1;
        ind(perm) = ind;
        
        testing_indices{ifold} = find(ind==0);
        training_indices{ifold} = find(ind==1);
        
    end
end

function [testing_indices, training_indices] = get_indices(ifold, k, n, direction)

    num = zeros(k, 1);  % number of entries in fold
    quotient = floor(n / k);
    remainder = mod(n, k);
    for i = 1:k
        num(i) = quotient;
        if (direction == 0 && i <= remainder) || (direction == 1 && i > k - remainder)
            num(i) = num(i) + 1;
        end
    end
    
    cumnum = cumsum(num);
    if ifold == 1
        first = 1;
    else
        first = cumnum(ifold-1)+1;
    end
    last = cumnum(ifold);
  
    testing_indices = first:last;
    training_indices = setdiff(1:n, testing_indices);
    
end