clear
close all

% uncomment one of the models below:
% g = @(m) (m-2).*(m+2).*(m+2)/8; m_true=2; label = 'tworoot';% nonlinear model with two roots
g = @(m) (m-2).*(m-2).*(m-2)/8; m_true=2; label='oneroot'; % nonlinear model, one root
%  g = @(m) tanh(m); m_true = 0; label='tanh';
g = @(m) 3*m-2; m_true = 1; label='linear';% linear model, one root
% g = @(m) 1; m_true=0; label='constant';

figure();
mplot = linspace(-10,10,500);
plot(mplot,g(mplot));
%%
% m_true = 1;
g_obs = g(m_true);
noise_coefficient = 0.1;

% estimate m and a noise coefficient.
m_accepted = 1.0;

steps = 5e6;
save_start = steps/2;
save_skip = 10;
isave = 1;
nsave = ceil((steps-save_start)/save_skip);
m_store = zeros(nsave,1);
g_store = zeros(nsave,1);
likeprob_store = zeros(nsave,1);

alpha_store = zeros(steps,1);
accept = false(steps,1);


iter=1;
while iter <= steps
    
    success = false;
    while ~success
        m_proposed = m_accepted;
        
        % perturb m
        m_proposed = m_proposed + randn()*0.1;
        % alternatively, pick m randomly across entire allowable range:
%         m_proposed = -10 + 20*rand();
        
        if m_proposed < -10 || m_proposed > 10
            success = false;
        else
            success = true;
        end
    end
    misfit_proposed = g(m_proposed)-g_obs;
    misfit_accepted = g(m_accepted)-g_obs;
    % likelihood of proposed and accepted solutions:
    %     likeprob_proposed = log(1) - 0.5*log(noise_coefficient*2*pi) + -0.5*misfit_proposed^2/noise_coefficient;
    likeprob_proposed = -0.5*misfit_proposed^2/noise_coefficient;
    %     likeprob_accepted = log(1) - 0.5*log(noise_coefficient*2*pi) + -0.5*misfit_accepted^2/noise_coefficient;
    likeprob_accepted = -0.5*misfit_accepted^2/noise_coefficient;
    % calculate acceptance probability - note that N=1
    
    alpha = likeprob_proposed - likeprob_accepted;
    alpha_store(iter) = alpha;
    if alpha > log(1) || log(rand()) <= alpha % accept/reject solution
        % accept solution
        m_accepted = m_proposed;
        likeprob_accepted = likeprob_proposed;
        accept(iter) = true;
    end
    if iter >= save_start && ~mod(iter,save_skip) % save results
        m_store(isave) = m_accepted;
        g_store(isave) = g(m_accepted);
        misfit_store(isave) = misfit_accepted;
        likeprob_store(isave) = likeprob_accepted;
        
        isave = isave + 1;
    end
    
    iter = iter + 1;
end

%% plotting
f=figure();
f.Position(3:4) = [560 600];
t = tiledlayout(4,3,'TileSpacing','none');
% t.Units = 'centimeters';
% % t.OuterPosition(3:4) = [16 16];

% nexttile(5,[2 2]);
[N,m_edges,g_edges] = histcounts2(m_store,g_store,linspace(-5,5,200),linspace(-5,5,201),'Normalization','pdf');
m_center = (m_edges(1:end-1)+m_edges(2:end))/2;
g_center = (g_edges(1:end-1)+g_edges(2:end))/2;
% pcolor(m_center,g_center,N');
% hcb = colorbar(); shading flat;
% hcb.Label.String = 'Probability Density';
% colormap(flipud(gray));
% % set(gca,'ColorScale','log');
% set(gca,'FontSize',12,'FontName','Helvetica');
% set(gca,'Box','on');
% text(0.05,0.90,'C','Units','normalized','FontSize',16,'FontName','Helvetica');
% ax2 = gca();
% set(gca,'layer','top','Box','on');
nexttile(5,[2 2]);
nplot=1e4;
indplot = randperm(length(m_store),nplot);
scatter(m_store(indplot),g_store(indplot),[],likeprob_store(indplot),'.');
hcb=colorbar();
hcb.Label.String = 'p(g(m)|d)'
set(gca,'FontSize',12,'FontName','Helvetica');
set(gca,'Box','on');
text(0.05,0.90,'C','Units','normalized','FontSize',16,'FontName','Helvetica');
ax2 = gca();

m_mean = mean(m_store); m_mean_style = 'r-';
r_mean = sqrt( (g(m_mean) - g_obs).^2 );
m_median = median(m_store); m_median_style = 'g--';
r_median = sqrt( (g(m_median) - g_obs).^2 );
r_true = sqrt( (g(m_true) - g_obs).^2 ); m_true_style = 'k-';

% Plot of m
nexttile(8+3,[1 2]);
histogram(m_store,m_edges,'EdgeColor','none','FaceColor',0.65*[1 1 1],'Normalization','pdf');
hold on
plot(m_mean*[1 1],get(gca,'YLim'),m_mean_style);
plot(m_median*[1 1],get(gca,'YLim'),m_median_style);
plot(m_true*[1 1],get(gca,'YLim'),m_true_style);

xlabel('m');
ylabel('p(m|d)');
set(gca,'FontSize',12,'FontName','Helvetica');
text(0.05,0.80,'D','Units','normalized','FontSize',16,'FontName','Helvetica');
ax1 = gca();

% plot of g(m)
nexttile(1+3,[2 1]);
% [n,c] = histcounts( (g_store-g_obs).^2/noise_coefficient );
% barh(n,c);
histogram(g_store,g_edges,'Orientation','horizontal','FaceColor',0.65*[1 1 1],'EdgeColor','none','Normalization','pdf');
ylabel('g(m)')
hold on
% plot(g_obs.^2*[1 1],get(gca,'YLim'),'k--');
plot(get(gca,'XLim'),g(m_mean)*[1 1],m_mean_style);
plot(get(gca,'XLim'),g(m_median)*[1 1],m_median_style);
plot(get(gca,'XLim'),g(m_true)*[1 1],m_true_style);
set(gca,'FontSize',12,'FontName','Helvetica');
text(0.1,0.9,'B','Units','normalized','FontSize',16,'FontName','Helvetica');
% xlabel('p(g(m)|d)');
ax3 = gca();
set(ax3,'YTick',get(ax2,'YTick'));
set(ax3,'YLim', get(ax2,'YLim'));
set(ax1,'XTick',get(ax2,'XTick'));
set(ax1,'XLim', get(ax2,'XLim'));
% set(ax2,'XTickLabel',[]);
% set(ax2,'YTickLabel',[]);
set(gcf,'Color','w');

% cumulative distribution
nexttile(1,[1 2]);
histogram(sqrt(misfit_store.^2),'Normalization','cdf','EdgeColor','none','FaceColor',0.65*[1 1 1]);
hold on
plot(r_true*[1 1],get(gca,'YLim'),m_true_style);
plot(r_mean*[1 1],get(gca,'YLim'),m_mean_style);
plot(r_median*[1 1],get(gca,'YLim'),m_median_style);
xlabel('||d-g(m)||_2');%,'interpreter','latex');
ylabel('cum. prob.');
set(gca,'FontSize',12,'FontName','helvetica');
text(0.05,0.85,'A','Units','normalized','FontSize',16,'FontName','Helvetica');

%% saving
save = 0
if save
    fig=gcf();
    set(fig,'Visible','off');
    set(fig,'renderer','painters');
    exportgraphics(t,[label '-figure.pdf']);
    set(fig,'renderer','opengl');
    set(fig,'Visible','on');
end
%% plot distribution of misfit
figure();
histogram(misfit_store);
hold on
% plot(r_mean*[1 1],get(gca,'YLim'),m_mean_style);
% plot(r_median*[1 1],get(gca,'YLim'),m_median_style);
%% calculate acceptance rate
N = 1000;
acceptance_rate = zeros(steps-N,1);
for i=1:(steps-N)
    acceptance_rate(i) = mean(accept(i:i+N-1));
end
figure, plot(acceptance_rate);
hold on
plot(get(gca,'XLim'),0.3*[1 1]);
