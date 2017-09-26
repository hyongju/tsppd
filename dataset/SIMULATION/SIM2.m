function [ output_args ] = SIM2( case_number )

n = length(case_number);
for i = 1:n,
   sample_num = case_number(i)
   simulation([10,20,40,50],sample_num); 
end

end

