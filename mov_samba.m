% Step 2: MOV calculation
%
% Load the Samba-W and Samba-E concatenated on Step 1

%load concat_IESsamba

%gpan_A= gsw_geo_strf_dyn_height(sal_A,tem_A,pres,1000);
gpan_AA= gsw_geo_strf_dyn_height(sal_AA,tem_AA,prew,1000);
gpan_B= gsw_geo_strf_dyn_height(sal_B,tem_B,prew,1000);
gpan_BB= gsw_geo_strf_dyn_height(sal_BB,tem_BB,prew,1000);
gpan_C= gsw_geo_strf_dyn_height(sal_C,tem_C,prew,1000);
gpan_CC= gsw_geo_strf_dyn_height(sal_CC,tem_CC,prew,1000);
gpan_D= gsw_geo_strf_dyn_height(sal_D,tem_D,prew,1000);

gpan_P1= gsw_geo_strf_dyn_height(sal_P1,tem_P1,pree,1000);
gpan_P2= gsw_geo_strf_dyn_height(sal_P2,tem_P2,pree,1000);
gpan_P4= gsw_geo_strf_dyn_height(sal_P4,tem_P4,pree,1000);
gpan_P5= gsw_geo_strf_dyn_height(sal_P5,tem_P5,pree,1000);
gpan_P6= gsw_geo_strf_dyn_height(sal_P6,tem_P6,pree,1000);
gpan_P8= gsw_geo_strf_dyn_height(sal_P8,tem_P8,pree,1000);

save gpan_samba gpan* lone lonw dt pre* depe*

%for ii=1:length(jd)
%ga=[gpanE(:,ii) gpanF(:,ii)];	 
%%vel(:,ii)=sw_gvel(ga,lat,lon);
%vel(:,ii)=gsw_geostrophic_velocity(ga,lon,lat,pres(:,1));
%end

%dist=gsw_distance(lo,lo,'km'); % distance im [m]

%%vref=repmat(vel(51,:),165,1);
%%vv=vel-vref;
%vv=vel(1:51,:); % ref @ 1000dbar
%transp=flipud(cumsum(flipud(vv)));
%plot(jd,transp(1,:)*20*dist/1e6)
%datetick('x',10);grid

%%title('Geostrophic transport between E and F, ref to 1000dbar')
%xlabel('year')
%ylabel('transport (Sv)')
%print -dpng transp_EF

%save transp_cpies
