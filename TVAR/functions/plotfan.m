function [ output_args ] = plotfan( fore )

% Set relevant quantiles (in percent)
quantiles = 5:5:95;
nq        = size(quantiles,2);

% Compute quantiles
mat_quant =  quantile(fore',quantiles/100)';

% Prepare plot matrix for its use with the area function
matm      = [mat_quant(:,1) mat_quant(:,2:end)-mat_quant(:,1:end-1)];


% Generate plot
h  = area(matm);
c  = .9:-.1:.1;
ct = [c sort(c)];

set(h,'LineStyle','none')
set(h(1), 'FaceColor', [1 1 1]) % white
for i = 1:nq-1
    set(h(i+1), 'FaceColor',[ct(i) ct(i)  ct(i)  ])
end
hold on
hist=median(fore');
plot(hist,'-r');  % median forecast in black

% set(gcf,'Color','w')


end

