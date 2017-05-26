% Script to generate Figure 4: comparing CACP method to CP method for the 
% unit sphere in R 3

% Select numbers of grid cells
M = [40,60,80,120,160,240,320];

% Specify number of threads to build system with
numthreads = 1;

% Do not generate results if already present
N = length(M);
if (exist('sphere_results.mat','file') ~=2 )
    CP_L2_error = zeros(1,N);
    CACP_L2_error = zeros(1,N);
    CP_LI_error = zeros(1,N);
    CACP_LI_error = zeros(1,N);
    CP_nnz = zeros(1,N);
    CACP_nnz = zeros(1,N);
    CP_condition = zeros(1,N);
    CACP_condition = zeros(1,N);
    for j = 1:N
        if (numthreads == 1)
            [e2,eI,numnz,condn,X,Y,Ucp,Ucacp,Uex] = shiftp_equation3(M(j));
        else
            [e2,eI,numnz,condn,X,Y,Ucp,Ucacp,Uex] = shiftp_equation3(M(j), numthreads);
        end
        CP_L2_error(j) = e2(1);
        CACP_L2_error(j) = e2(2);
        CP_LI_error(j) = eI(1);
        CACP_LI_error(j) = eI(2);
        CP_nnz(j) = numnz(1);
        CACP_nnz(j) = numnz(2);
        CP_condition(j) = condn(1);
        CACP_condition(j) = condn(2);
    end
    
    save('sphere_results.mat','CP_L2_error','CACP_L2_error','CP_LI_error',...
        'CACP_LI_error','CP_nnz','CACP_nnz','CP_condition','CACP_condition');
else
    load('sphere_results.mat');
end   

lw = 2;     %   linewidth
ms = 10;    %   marker size
fs = 16;    %   font size

figure(1)
clf
plot(M,CACP_L2_error./CP_L2_error,'-bo',...
     M,CACP_LI_error./CP_LI_error,'-rs','linewidth',lw,'markersize',ms);
xlabel('number of grid cells (M)','Fontweight','bold','fontsize',fs);
ylabel('normalized error','Fontweight','bold','fontsize',fs);
set(gca,'fontsize',fs);
outerpos = get(gca,'OuterPosition');
ti = get(gca,'TightInset');
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
set(gca,'Position',[left bottom ax_width, ax_height]);
saveas(gcf,'normalized_error_sphere','png');

figure(2)
clf
n = length(M);
aCP = sum( (CP_nnz(2:end) - CP_nnz(1:end-1))./(M(2:end).^2 - M(1:end-1).^2) )/(N-1);
bCP = sum(CP_nnz - aCP*M.^2)/N;
aCACP = sum( (CACP_nnz(2:end) - CACP_nnz(1:end-1))./(M(2:end).^2 - M(1:end-1).^2) )/(N-1);
bCACP = sum(CACP_nnz - aCACP*M.^2)/N;
plot(M,CP_nnz,'bo',M,aCP*M.^2 + bCP,'-b',M,CACP_nnz,'rs',M,aCACP*M.^2 + bCACP,'-r','linewidth',lw,'markersize',ms);
xlabel('number of grid cells (M)','Fontweight','bold','fontsize',fs);
ylabel('number of non-zero matrix entries','Fontweight','bold','fontsize',fs);
set(gca,'fontsize',fs);
outerpos = get(gca,'OuterPosition');
ti = get(gca,'TightInset');
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
set(gca,'Position',[left bottom ax_width, ax_height]);
saveas(gcf,'num_nonzeros_sphere','png');
