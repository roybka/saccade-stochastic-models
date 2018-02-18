function [ chi_square ] = test_fit2( orig,model )
% takes gaps of two distributions, and returns chisquare
if isempty(model)|isempty(orig)
    chi_square=inf;
    return
end
original_distribution=orig(orig<1000);
amount_of_bins=round(1.88*power(length(original_distribution),2/5));

%finding the numbers that devide it into equiprobable quantiles
bin_borders=[quantile(original_distribution,amount_of_bins-1)];

bin_borders=[0,bin_borders,1000]; %adding the borders of the bins;
%histc, should give the info of how many samples fall inside each bin
bincounts_original=histc(original_distribution,bin_borders)';
%edit it so the info will be the same for x and bincounts_original
bincounts_original=bincounts_original(1:end-1);
original_bin_percent=bincounts_original/sum(bincounts_original);

test_distribution=model;


%trim the simulated train between 0 and 1000;
test_distribution=test_distribution(test_distribution>0 & test_distribution<1000);

%using the borders from the original samples, and counting the bin counts
%in the simulated one:
bincounts_simulated=histc(test_distribution,bin_borders)';
bincounts_simulated=bincounts_simulated(1:end-1);
simulated_bin_percent=bincounts_simulated/sum(bincounts_simulated);

chi_square=sum(power((simulated_bin_percent-original_bin_percent),2)./original_bin_percent);

end

