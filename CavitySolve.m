% CavitySolve computes panel solution for constant strength source and
% dipoles distributed on panels of body.

clear

U = 1;
Udir = 0 * (pi/180);
Qinf = U * [cos(Udir) sin(Udir)];
alphabody = 30;

c = 1;          % chord
lcav = 3*c;     % cavity length;
Panels = 200;   % # of Panels on each sub-body
deres = 0;      % cavity panel spacing from detachment to closure, 1 = half cosine (low res at closure), 0 = full cosine (high res at closure)
cavshape = 3;   % cavity shape: 0 = linear, 1 = two single quadrant ellipse, 2 = 4-quadrant ellipse (Careful), 3 = 3-quadrant ellipse
cavoff = -.0;   % fractional vertical offset of cavity closure relative to width for initial cavity shape {-1 = 100% shift down, +1 = 100% shift up}
relax = 1;      % relaxation factor for cavity correction
clamped = 0;    % if clamped {=1}, cavity correction is zero at closure point
perturb = 1;    % if perturbed {=1}, perturb vertical cavity location to stabilize closure
its = 4;

thN1 = [1,0];
Lwake = 20*lcav;

% [wetlo,wethi,cavlo,cavhi] = Body_Cav_Composite(Panels,alphabody,c,lcav,deres,cavshape,cavoff);
% wet = [flipud(wetlo(2:end,:));wethi];
% Lrefhi = Pan(end,1)-wethi(end,1);
% Lreflo = Pan(1,1)-wetlo(end,1);

wet = Bodyshape(Panels,alphabody,c,2,1,1);
[cavlo,cavhi] = Cavityshape(wet,lcav,c);
cavlo = flipud(cavlo);
Pan = [cavlo(1:end-1,:);wet(1:end-1,:);cavhi];
Lrefhi = Pan(end,1)-wet(end,1);
Lreflo = Pan(1,1)-wet(1,1);

Nw = length(wet)-1;
Ncl = length(cavlo)-1;
Nch = length(cavhi)-1;
N = Nw + Ncl + Nch;

Swet = lengths(wet);
% SdLo = CubExtrap(Swet,0);
% SdHi = fliplr(CubExtrap(flipud(Swet),0));
SdLo = Extrap(Swet,0,4);
SdHi = fliplr(Extrap(flipud(Swet),0,4));

sig = zeros(N,its);                      % interior potential referenced INTO body
mu = zeros(N,its);
qc = zeros(1,its);

skew = -0.1:0.1:0.1;
fixset = zeros(1,length(skew));

for now = 1:its

    PanTemp = Pan;
    PanOld = Pan;
    
    if perturb == 1

        for j = 1:length(skew)

            PanTemp(1:Ncl+1,:)=skewy(Pan(1:Ncl+1,:),skew(j),1);
            PanTemp(N+1-Ncl:end,:)=skewy(Pan(N+1-Ncl:end,:),skew(j),0);

            [nh,th] = nhat(PanTemp,1);
            xN1 = PanTemp(1,:) + Lwake*thN1;
            Scl = lengths(PanTemp(1:Ncl+1,:));
            Sch = lengths(PanTemp(Ncl+1+Nw:N+1,:));

            [Big,Rhs] = BuildSysPar(PanTemp,N,Ncl,Nch,Nw,SdLo,SdHi,Scl,Sch,xN1,Qinf,clamped,nh);

            soln = Big\Rhs;
            mu(:,now) = soln(1:N);
            sig(:,now) = soln(N+1:2*N);
            qc(now) = soln(2*N+1);

            [Pan2,fixset(j)] = Corrections(PanTemp,nh,sig(:,now),qc(now),Ncl,Nch,Nw,Scl,Sch,Qinf,relax);

%             figure;
%             hold on
%             plot(Pan(:,1),Pan(:,2),'bl');
%             plot(Pan(:,1),Pan(:,2),'bl.');
%             plot(Pan2(:,1),Pan2(:,2),'k.');
%             plot(Pan2(:,1),Pan2(:,2),'g');
%             plot(Pan2(1,1),Pan2(1,2),'ro');
%             plot(Pan2(Ncl,1),Pan2(Ncl,2),'bo');
%             plot(Pan2(Ncl+Nw+2,1),Pan2(Ncl+Nw+2,2),'ko');
%             plot(Pan2(N+1,1),Pan2(N+1,2),'go');
%             axis equal;
%             grid on;
            
        end;

        ply = polyfit(skew,fixset,2);
        g = @(x) polyval(ply,x);
        shootfor = fzero(g,0);

        Pan(1:Ncl+1,:)=skewy(Pan(1:Ncl+1,:),shootfor,1);
        Pan(N+1-Ncl:end,:)=skewy(Pan(N+1-Ncl:end,:),shootfor,0);
        
    end;

    [nh,th] = nhat(Pan,1);
    xN1 = Pan(1,:) + Lwake*thN1;
    Scl = lengths(Pan(1:Ncl+1,:));
    Sch = lengths(Pan(Ncl+1+Nw:N+1,:));
    
    PanOld = Pan;
    xN1Old = xN1;

    [Big,Rhs] = BuildSysPar(Pan,N,Ncl,Nch,Nw,SdLo,SdHi,Scl,Sch,xN1,Qinf,clamped,nh);

    soln = Big\Rhs;
    mu(:,now) = soln(1:N);
    sig(:,now) = soln(N+1:2*N);
    qc(now) = soln(2*N+1);
    

    [Pan2,blank,hch,hcl] = Corrections(Pan,nh,sig(:,now),qc(now),Ncl,Nch,Nw,Scl,Sch,Qinf,relax);
    
    Pan = FixPanCubic(Pan,Pan2,Lrefhi,Lreflo,N,Ncl,Nch,Nw);
%     Pan = FixPanArc(Pan,Pan2,Lrefhi,Lreflo,N,Ncl,Nch,Nw);

    [nh,th] = nhat(Pan,1);
    xN1 = Pan(1,:) + Lwake*thN1;
    Scl = lengths(Pan(1:Ncl+1,:));
    Sch = lengths(Pan(Ncl+1+Nw:N+1,:));
        
%     figure;
%     hold on
%     plot(Pan(:,1),Pan(:,2),'k');
%     quiver(Pan(1:N,1),Pan(1:N,2),nh(:,1),nh(:,2),'r');
%     quiver(Pan(1:N,1),Pan(1:N,2),th(:,1),th(:,2),'g');
%     plot([Pan(N+1,1),xN1(1)],[Pan(N+1,2),xN1(2)],'r');
%     plot(cavlo(:,1),cavlo(:,2),'c');
%     plot(cavhi(:,1),cavhi(:,2),'m');
%     axis equal;
%     grid on;
%     xlim([-c,1.5*lcav]);
    
%     figure;
%     hold on;
%     plot(flipud(sig(1:Ncl,now)),'g.');
%     plot(sig(Ncl+Nw+1:end,now),'r.');
%     plot(hcl,'g');
%     plot(hch,'r');
%     grid on;

end;

% figure;
% hold on;
% plot(sig,'r');
% plot(mu,'k');

%%
[Big,Rhs] = BuildSysPar(Pan,N,Ncl,Nch,Nw,SdLo,SdHi,Scl,Sch,xN1,Qinf,clamped,nh);

PanOld = Pan;
xN1Old = xN1;

soln = Big\Rhs;
mu(:,now) = soln(1:N);
sig(:,now) = soln(N+1:2*N);
qc(now) = soln(2*N+1);

% [Pan] = Corrections(Pan,nh,sig(:,now),qc(now),Ncl,Nch,Nw,Scl,Sch,Qinf,relax);
% [nh,th] = nhat(Pan,1);
% xN1 = Pan(1,:) + 100*lcav*thN1;

% sect=Pan(2:end,:)-Pan(1:end-1,:);
span=lengths(Pan);
span=(span(1:end-1)+span(2:end))/2;
Qtinf=Qinf*th'; 
Qt=-(mu(2:end,end)-mu(1:end-1,end))./span + Qtinf(1:end-1)';
figure;
plot(Qt,'k.')


%%
mu = [mu;mu(end,:)-mu(1,:)];
% mids = (Pan(1:N,:)+Pan(2:N+1,:))/2;
% midwake = (xN1+Pan(N+1,:))/2;
% mids = [mids;midwake];
% sects = (Pan(2:N+1,:)-Pan(1:N,:));
% sectwake = xN1-Pan(N+1,:);
% sects = [sects;sectwake];
% lengths = sqrt(sects(:,1).^2+sects(:,2).^2);
% sigd = [lengths(1:N).*sig(:,end);0];
% mud = lengths.*mu(:,end);
% alpha_n = [atan2(nh(:,2),nh(:,1));atan2(xN1(1)-Pan(N+1,1),Pan(N+1,2)-xN1(2))];
% gamma = mu(2:N+1) - mu(1:N);
% gamma = [gamma; -gamma(N)];

%%
figure;
plot(Pan(1:N+1,1),mu(:,end),'k.');
hold on
plot(Pan(N+1,1),mu(N+1,end),'ro');
plot(Pan(1:N,1),sig(:,end),'g.');


% Qt = -(mu(2:N+1,end)-mu(1:N,end))./lengths(1:N) + (th*Qinf');
Cp = 1-Qt.^2./(Qinf*Qinf');
figure;
plot(Pan(2:N,1),-Cp(1:N-1),'k.');
hold on;
plot(Pan(1:N,1),Pan(1:N,2));
% figure;
% plot(Pan(2:N,1),Qt(1:N-1));


[X,Y,Z,Zp,Zm]=StreamPlots([-1.5*c lcav+c ; -c c],0.05,Qinf,PanOld,sig(:,end),mu(:,end));

[X,Y,Z,Zp,Zm]=StreamPlots([-1.5*lcav 2*lcav ; -lcav lcav],0.05,Qinf,PanOld,sig(:,end),mu(:,end));



