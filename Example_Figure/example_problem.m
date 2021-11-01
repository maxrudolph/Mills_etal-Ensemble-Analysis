clear
close all

% uncomment one of the models below:
% g1 = @(m) (m-2).*(m+2).*(m+2)/8; m_true=-3; label = 'tworoot';% nonlinear model with two roots
g1 = @(m) (m-2).*(m-2).*(m-2)/8; m_true=1; label='oneroot'; % nonlinear model, one root
%  g1 = @(m) tanh(m); m_true = 0; label='tanh';
% g1 = @(m) m; m_true = 1; label='linear';% linear model, one root
% g1 = @(m) 1; m_true=0; label='constant';

% figure();
% mplot = linspace(-10,10,500);
% plot(mplot,g(mplot));

% introduce an outer loop over m-values
all_m = [1 -0.5];
plotcolors = colororder;

for im = 1:length(all_m)
    %% initialize parameters
    m_true = all_m(im);
    
    g = @(m) g1(m);% + 2*(im-1); % add an offset for plotting beauty
    g_obs = g(m_true);
    
    noise_coefficient = 0.1; % this is sigma^2
    
    m_accepted = 1.0;
    
    steps = 5e8;
    save_start = steps/2;
    save_skip = 100;
    isave = 1;
    nsave = ceil((steps-save_start)/save_skip);
    m_store = zeros(nsave,1);
    g_store = zeros(nsave,1);
    likeprob_store = zeros(nsave,1);
    misfit_store = zeros(nsave,1);
    
    alpha_store = zeros(steps,1);
    accept = false(steps,1);
    
    likeprob_accepted = -1e99;
    iter=1;
    hwb = waitbar(0,'Progress');
    while iter <= steps
        
        success = false;
        while ~success
            m_proposed = m_accepted;
            
            % perturb m
            m_proposed = m_proposed + randn()*0.05;
            % alternatively, pick m randomly across entire allowable range:
            %         m_proposed = -10 + 20*rand();
            
            if m_proposed < -10 || m_proposed > 10
                success = false;
            else
                success = true;
            end
        end
        misfit_proposed = g(m_proposed)-g_obs;
        %     misfit_accepted = g(m_accepted)-g_obs;
        % likelihood of proposed and accepted solutions:
        %     likeprob_proposed = log(1) - 0.5*log(noise_coefficient*2*pi) + -0.5*misfit_proposed^2/noise_coefficient;
        likeprob_proposed = -0.5*misfit_proposed^2/noise_coefficient;
        %     likeprob_accepted = log(1) - 0.5*log(noise_coefficient*2*pi) + -0.5*misfit_accepted^2/noise_coefficient;
        %     likeprob_accepted = -0.5*misfit_accepted^2/noise_coefficient;
        % calculate acceptance probability - note that N=1
        
        alpha = likeprob_proposed - likeprob_accepted;
        alpha_store(iter) = alpha;
        if alpha > log(1) || log(rand()) < alpha % accept/reject solution
            % accept solution
            m_accepted = m_proposed;
            misfit_accepted = misfit_proposed;
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
        if ~mod(iter,10000)
            waitbar(iter/steps,hwb);
        end
        iter = iter + 1;
    end
    close(hwb);
    
    %% analytic calculation
    % p(m|I)
    prior = @(m) 1/20; % 1/(10- -10)
    % p(d|m,I)
    prob = @(m) 1/sqrt(noise_coefficient*2*pi) * exp(-0.5*(g(m)-g_obs).^2/noise_coefficient);
    evidence = @(m) 1.;% p(d|I)
    
    mvals = linspace(-5,5,1e6);
    posterior = prior(mvals).*prob(mvals)./evidence(mvals); % this is const*p(m|d,I).
    cdf = cumtrapz(mvals,posterior);
    
    % remove duplicate values of 0/1 at end of cdf
    mask = find(cdf == 0,1,'last'):find(cdf == cdf(end),1,'first');
    normalization_factor = cdf(end);
    posterior = posterior/normalization_factor;
    
    nsamples = 1e7;
    tmp = rand(nsamples,1);
    m_samples = interp1(cdf(mask),mvals(mask),tmp);
    
    figure();
    subplot(2,1,1);
    % histogram(m_samples,'EdgeColor','none','Normalization','pdf');
    [f,xi] = ksdensity(m_samples);
    plot(xi,f,'b');
    hold on
    plot(mvals,posterior,'r--');
    
    xlabel('m');
    subplot(2,1,2);
    histogram(g(m_samples),'Normalization','pdf');
    xlabel('g');       
    
    %% plotting
    if im == 1
        f=figure(901);
        clf();
        f.Position(3:4) = [560 600];
        t = tiledlayout(4,3,'TileSpacing','none');
    else
        f=figure(901);
    end
    alpha = 0.15;
    
    
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
    nplot=5e3;
    indplot = randperm(length(m_store),nplot);
    
    if im==1
        mtmp = linspace(-1,4,201);
        plot(mtmp,g(mtmp),'Color','k');%[plotcolors(im,:) 0.25]);
        hold on
    end
    colors = [plotcolors(im,:); 1 1 1];
    colors = interp1([min(-likeprob_store) max(-likeprob_store)],colors,-likeprob_store(indplot));
    sizes = interp1([min(-likeprob_store) max(-likeprob_store)],[20 0],-likeprob_store(indplot));
    hs=scatter(m_store(indplot),g_store(indplot),sizes,colors,'filled',...
         'Marker','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.01);       
    hold on;
    
    set(gca,'FontSize',12,'FontName','Helvetica');
    set(gca,'Box','on');
    text(0.05,0.90,'C','Units','normalized','FontSize',16,'FontName','Helvetica');
    ax2 = gca();
    
    lw=1; % line width
    m_mean = mean(m_store); m_mean_style = ':';
    r_mean = sqrt( (g(m_mean) - g_obs).^2 );
    m_median = median(m_store); m_median_style = '--';
    r_median = sqrt( (g(m_median) - g_obs).^2 );
    r_true = sqrt( (g(m_true) - g_obs).^2 ); m_true_style = '-';
    
    %
    % Panel D - Plot of m
    % 
    ymax = 3.2;
    nexttile(8+3,[1 2]);
    % [f,xi] = ksdensity(m_store);
    histogram(m_store,m_edges,'EdgeColor','none','FaceColor',plotcolors(im,:),'Normalization','pdf','FaceAlpha',alpha);
    % plot(xi,f);
    hold on
    [f,xi] = ksdensity(m_samples);
    plot(xi,f,'Color',plotcolors(im,:));
    plot(m_mean*[1 1],[0 ymax],m_mean_style,'Color',plotcolors(im,:),'LineWidth',lw);
    plot(m_median*[1 1],[0 ymax],m_median_style,'Color',plotcolors(im,:),'LineWidth',lw);
    plot(m_true*[1 1],[0 ymax],m_true_style,'Color',plotcolors(im,:),'LineWidth',lw);
    set(gca,'YLim',[0 ymax]);
    xlabel('m');
    ylabel('p(m|d)');
    set(gca,'FontSize',12,'FontName','Helvetica');
    text(0.05,0.80,'D','Units','normalized','FontSize',16,'FontName','Helvetica');
    ax1 = gca();
    
    %
    % Pabel B - plot of g(m)
    %
    nexttile(1+3,[2 1]);
    % [n,c] = histcounts( (g_store-g_obs).^2/noise_coefficient );
    % barh(n,c);
    histogram(g_store,g_edges,'Orientation','horizontal','FaceColor',plotcolors(im,:),'EdgeColor','none','Normalization','pdf','FaceAlpha',alpha);
    %     [f,xi] = ksdensity(g_store);
    %     plot(f,xi);
    hold on
    % histogram(g(m_samples),'Orientation','horizontal','FaceColor',0.65*[0 0 1],'EdgeColor','none','Normalization','pdf');
    [f,xi] = ksdensity(g(m_samples));
    plot(f,xi,'Color',plotcolors(im,:),'LineWidth',lw);
    ylabel('g(m)')
    hold on
    % plot(g_obs.^2*[1 1],get(gca,'YLim'),'k--');
    plot(get(gca,'XLim'),g(m_mean)*[1 1],m_mean_style,'Color',plotcolors(im,:),'LineWidth',lw);
    plot(get(gca,'XLim'),g(m_median)*[1 1],m_median_style,'Color',plotcolors(im,:),'LineWidth',lw);
    plot(get(gca,'XLim'),g(m_true)*[1 1],m_true_style,'Color',plotcolors(im,:),'LineWidth',lw+1);
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
    nexttile(1,[1 3]);
    histogram(sqrt(misfit_store.^2),'Normalization','cdf','EdgeColor','none','FaceColor',plotcolors(im,:),'FaceAlpha',alpha);
    hold on    
    plot(r_true*[1 1],get(gca,'YLim'),m_true_style,'Color',plotcolors(im,:),'LineWidth',lw+1);
    plot(r_mean*[1 1],get(gca,'YLim'),m_mean_style,'Color',plotcolors(im,:),'LineWidth',lw);
    plot(r_median*[1 1],get(gca,'YLim'),m_median_style,'Color',plotcolors(im,:),'LineWidth',lw);
    [f,xi] = ksdensity( sqrt((g(m_samples)-g_obs).^2 ),'Function','cdf');
    plot(xi,f,'Color',plotcolors(im,:));
    xlabel('||d-g(m)||_2');%,'interpreter','latex');
    ylabel('cum. prob.');
    set(gca,'FontSize',12,'FontName','helvetica');
    set(gca,'XLim',[-.05 1]);
    text(0.05,0.85,'A','Units','normalized','FontSize',16,'FontName','Helvetica');
    
    %% saving    
    save = 1
    if save && im == length(all_m)
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
    % calculate acceptance rate
    
    acceptance_rate = movmean(accept,1000);
    figure, plot(acceptance_rate);
    hold on
    plot(get(gca,'XLim'),0.3*[1 1]);
    
end
