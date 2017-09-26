function [ output_args ] = SIM3( case_number )

n = length(case_number);
for i = 1:n,
   sample_num = case_number(i)
   simulation([12,30],sample_num); 
end

end

