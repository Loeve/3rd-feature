function [distance] = chi_square_distance(vectorA, vectorB)

vectorA = double(vectorA);
vectorB = double(vectorB);

vectorA = vectorA ./ sum(vectorA);
vectorB = vectorB ./ sum(vectorB);

distance = 0.0;
squared_difference = (vectorA - vectorB).^2;
sum_vector = vectorA + vectorB;
lengthA = length(vectorA);

for i = 1: lengthA
    denominator = sum_vector(i);
    if (denominator > 0)
        distance = distance + (squared_difference(i) / denominator);
    end
end

distance = 0.5 * distance;

end

