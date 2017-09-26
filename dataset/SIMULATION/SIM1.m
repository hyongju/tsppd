function [ output_args ] = SIM1( case_number )

n = length(case_number);
for i = 1:n,
   sample_num = case_number(i)
   simulation([8,15,25],sample_num); 
end

end

