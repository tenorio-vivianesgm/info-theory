function H=JointEntropy_byJuli(matrix)

% matrix is the realization of a m (columns) dimensional random variable
% with sample size n (number of rows)

[n m] = size(matrix);

Alphabet = unique(matrix, 'rows');
Frequency = zeros(size(Alphabet,1),1);

for sample = 1:n
    for symbol=1:size(Alphabet, 1)
        if matrix(sample,:)== Alphabet(symbol,:)
            Frequency(symbol) = Frequency(symbol)+1;
        end
    end
end

P = Frequency / sum(Frequency);
H = -sum(P .* log2(P));