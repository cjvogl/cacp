% Script to generate Figure 2: CACP results for the clover in R^2

% Select numbers of grid cells
M = [40,60,80,120,160,240,320,480,640,960];

% Do not generate results if already present
N = length(M);
if (exist('clover_results.mat','file') ~=2 )
    CP_L2_error = zeros(1,N);
    CACP_L2_error = zeros(1,N);
    CP_LI_error = zeros(1,N);
    CACP_LI_error = zeros(1,N);
    CP_nnz = zeros(1,N);
    CACP_nnz = zeros(1,N);
    CP_condition = zeros(1,N);
    CACP_condition = zeros(1,N);
    for j = 1:N
        [e2,eI,numnz,condn,X,Y,Ucp,Ucacp,Uex] = shiftp_equation(M(j), 2);
        CP_L2_error(j) = e2(1);
        CACP_L2_error(j) = e2(2);
        CP_LI_error(j) = eI(1);
        CACP_LI_error(j) = eI(2);
        CP_nnz(j) = numnz(1);
        CACP_nnz(j) = numnz(2);
        CP_condition(j) = condn(1);
        CACP_condition(j) = condn(2);
    end
    
    save('clover_results.mat','CP_L2_error','CACP_L2_error','CP_LI_error',...
        'CACP_LI_error','CP_nnz','CACP_nnz','CP_condition','CACP_condition');
else
    load('clover_results.mat');
end    

lw = 2;     %   linewidth
ms = 10;    %   marker size
fs = 16;    %   font size

figure(1)
clf
[e2,eI,numnx,condn,x,y,uold,unew,uex] = shiftp_equation(160, 2);
unew(unew == 0) = nan;
pcolor(x,y,unew)
colorbar
colormap copper
shading interp
axis equal
set(gca,'fontsize',fs);
outerpos = get(gca,'OuterPosition');
ti = get(gca,'TightInset');
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3) - 0.15;
ax_height = outerpos(4) - ti(2) - ti(4);
set(gca,'Position',[left bottom ax_width, ax_height]);
saveas(gcf,'solution_clover','png');

figure(2)
clf
dx = 2./M(2:2:end);
loglog(dx,CACP_L2_error(2:2:end),'-bo',dx,CACP_LI_error(2:2:end),'-sr',dx,10*dx.^2,'k--','linewidth',lw,'markersize',ms)
xlabel('{\Delta}x','Fontweight','bold','fontsize',fs)
ylabel('error','Fontweight','bold','fontsize',fs)
set(gca,'fontsize',fs);
outerpos = get(gca,'OuterPosition');
ti = get(gca,'TightInset');
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
set(gca,'Position',[left bottom ax_width, ax_height]);
saveas(gcf,'convergence_clover','png');