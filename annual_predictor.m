function [prediction] = annual_predictor(previous_data,n_years)

num = 5;
% Categorising data
limits = [min(previous_data) quantile(previous_data,num) max(previous_data)+0.01];
y_cat = ones(1,length(previous_data));
for i = 1:length(previous_data)
    for j = 1:num+1
        if previous_data(i) >= limits(j) && previous_data(i) < limits(j+1)
            y_cat(i) = j;
        end
    end
end
        
% Transition matrix
transition_matrix = zeros(num+1,num+1);

for i = 1:length(previous_data)-1
    for j = 1:num+1
        for k = 1:num+1
            if y_cat(i) == j && y_cat(i+1) == k
                transition_matrix(j,k) = transition_matrix(j,k) + 1;
            end
        end
    end
            
end

sums = sum(transition_matrix,2);
trans_mat = zeros(num+1,num+1);

for i = 1:num+1
    trans_mat(i,:) = transition_matrix(i,:)./sums(i);
end

% Iteration
final_state = zeros(1,n_years);
prediction = zeros(1,n_years);

for i = 1:n_years
    l = randi(num+1);
    a = zeros(1,num+1);
    a(l) = 1;
    x = a * trans_mat;
    b = find(x == max(x));
    if length(b) > 1
        b = b(randi(length(b)));
    end
    final_state(i) = b;
    for j = 1:num+1
        if final_state(i) == j
            prediction(i) = limits(j) + (limits(j+1) - limits(j))*rand();
        end
    end
end
end

