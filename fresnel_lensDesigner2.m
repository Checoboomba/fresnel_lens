%Matlab Script to design a fresnel lens, veify it and a .stl file.
clear;clc;close all;
freq_GHz=valid_input('Enter operating frequency [220-330]: ',220,330);
diameter_mm=valid_input('Enter physical diameter in mm [100-200]: ',100,200);
active_diameter_mm=valid_input(sprintf('Enter optically active diameter in mm [must be < %.1f]: ',diameter_mm),1,diameter_mm-0.1);
f_work=valid_input('Enter focal length in mm: ',1,1000);
rim_width=(diameter_mm-active_diameter_mm)/2;
fprintf('\nMounting rim width: %.2f mm on each side\n',rim_width);
if rim_width<1.0
fprintf('Note! rim is very narrow (<1mm).\n');
end
n=1.521;
c_mm_s=3e11;
base_thickness=2.0;
lambda=c_mm_s/(freq_GHz*1e9);
rho_max=active_diameter_mm/2;
blank_radius=diameter_mm/2;
t_groove=lambda/(n-1);
d_ref=t_groove+base_thickness;
fprintf('\n--- Design summary ---\n');
fprintf('Wavelength: %.4f mm\n',lambda);
fprintf('Groove depth: %.4f mm\n',t_groove);
fprintf('Back-face position: %.4f mm\n',d_ref);
Npts=300;
rho=linspace(0,rho_max,Npts);
D=sqrt(f_work^2*(n-1)^2+(n^2-1).*rho.^2);
s_ideal=(-f_work*(n-1)+D)./(n^2-1);
m_zone=floor(s_ideal/t_groove);
s_fresnel=s_ideal-m_zone.*t_groove;
n_zones=max(m_zone);
fprintf('Number of Fresnel zones: %d\n',n_zones);
fprintf('Conventional-lens max sag: %.3f mm -> Fresnel max thickness: %.3f mm\n',max(s_ideal),max(s_fresnel));
s_prime=rho./D;
Nz=-1*ones(size(rho));
Nr=s_prime;
Nmag=sqrt(Nz.^2+Nr.^2);
Iz=s_fresnel+f_work;
Ir=rho;
Imag=sqrt(Iz.^2+Ir.^2);
theta1=acos(abs((Iz.*Nz+Ir.*Nr)./(Imag.*Nmag)));
theta2=acos(abs(Nz./Nmag));
snell_dev=abs(sin(theta1)-n.*sin(theta2));
d1=sqrt((s_fresnel+f_work).^2+rho.^2);
d2=d_ref-s_fresnel;
OPL=d1+n.*d2;
C0=f_work+n*d_ref;
implied_m=(C0-OPL)/lambda;
phase_err_deg=abs(implied_m-round(implied_m)).*360;
max_err=max(phase_err_deg);
theta_edge=atan(rho_max/(f_work+s_ideal(end)));
NA=sin(theta_edge);
fprintf('\n--- Verification ---\n');
fprintf('Numerical Aperture (NA): %.4f (f/# = %.3f)\n',NA,1/(2*NA));
fprintf('Max Snell deviation: %.4f\n',max(snell_dev));
fprintf('Max OPL phase error: %.1f degrees\n',max_err);
if max_err>90
fprintf('WARNING: exceeds Rayleigh quarter-wave criterion (90 deg).\n');
elseif max_err>45
fprintf('CAUTION: moderate phase error at outer zones.\n');
else
fprintf('OK: within acceptable performance.\n');
end
rho_full=[-fliplr(rho),rho];
s_full=[fliplr(s_fresnel),s_fresnel];

fig1=figure('Position',[100 100 850 750]);
hold on;grid on;axis equal;
plot(s_full,rho_full,'b-','LineWidth',1.3,'DisplayName','Active Fresnel zone');
plot([0 0],[rho_max blank_radius],'Color',[0.5 0.5 0.5],'LineWidth',2.5,'DisplayName','Flat mounting rim');
plot([0 0],[-blank_radius -rho_max],'Color',[0.5 0.5 0.5],'LineWidth',2.5,'HandleVisibility','off');
plot([0 d_ref],[blank_radius blank_radius],'--','Color',[0.5 0.5 0.5],'LineWidth',0.8,'HandleVisibility','off');
plot([0 d_ref],[-blank_radius -blank_radius],'--','Color',[0.5 0.5 0.5],'LineWidth',0.8,'HandleVisibility','off');
plot([d_ref d_ref],[-blank_radius blank_radius],'k-','LineWidth',2,'DisplayName','Flat exit face');
plot(-f_work,0,'rp','MarkerSize',15,'MarkerFaceColor','r','DisplayName','Source');
ray_heights=linspace(-rho_max*0.95,rho_max*0.95,11);
for rp=ray_heights
[~,idx]=min(abs(rho-abs(rp)));
sp=s_fresnel(idx);
plot([-f_work sp],[0 rp],'r-','LineWidth',0.6,'Color',[1 0.4 0.4],'HandleVisibility','off');
plot([sp d_ref],[rp rp],'-','LineWidth',1.0,'Color',[1 0.6 0],'HandleVisibility','off');
end
xlabel('z - optical axis (mm)');ylabel('\rho - aperture height (mm)');
title(sprintf('Fresnel Lens Focal Raytrace: %.0f GHz, f=%.0fmm',freq_GHz,f_work));
legend('Location','northwest','FontSize',8);
hold off;
exportgraphics(fig1,'fresnel_lens_raytrace.png','Resolution',150);

fig2=figure('Position',[150 150 950 650]);
hold on;grid on;axis equal;
plot(s_full,rho_full,'b-','LineWidth',1.3,'DisplayName','Fresnel Surface');
plot([0 0],[rho_max blank_radius],'Color',[0.5 0.5 0.5],'LineWidth',2.5,'DisplayName','Mounting Rim');
plot([0 0],[-blank_radius -rho_max],'Color',[0.5 0.5 0.5],'LineWidth',2.5,'HandleVisibility','off');
plot([d_ref d_ref],[-blank_radius blank_radius],'k-','LineWidth',2,'DisplayName','Flat Exit Face');
plot(-f_work,0,'rp','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Source Point');
z_exit_max=d_ref+f_work*1.5;
for rp=ray_heights
[~,idx]=min(abs(rho-abs(rp)));
sp=s_fresnel(idx);
plot([-f_work sp],[0 rp],'r-','LineWidth',0.6,'Color',[1 0.3 0.3],'HandleVisibility','off');
plot([sp d_ref],[rp rp],'-','LineWidth',1.0,'Color',[1 0.6 0],'HandleVisibility','off');
plot([d_ref z_exit_max],[rp rp],'g-','LineWidth',1.2,'Color',[0 0.7 0.2],'HandleVisibility','off');
end
plot([d_ref z_exit_max],[0 0],'k--','LineWidth',0.8,'DisplayName','Collimated Beam Axis');
xlabel('z - optical axis (mm)');ylabel('\rho - aperture height (mm)');
title(sprintf('Transmitted Ray Path Post-Lens (%.0f GHz Collimated Beam)',freq_GHz));
legend('Location','northwest','FontSize',8);
xlim([-f_work*1.1 z_exit_max*1.05]);
hold off;
exportgraphics(fig2,'fresnel_lens_exit_rays.png','Resolution',150);
fprintf('\nGenerated Ray-trace Images:\n 1) %s\n 2) %s\n',fullfile(pwd,'fresnel_lens_raytrace.png'),fullfile(pwd,'fresnel_lens_exit_rays.png'));

n_theta=100;
theta=linspace(0,2*pi,n_theta+1);theta(end)=[];
interior_rho=[rho(2:end),blank_radius,blank_radius];
interior_z=[s_fresnel(2:end),0,d_ref];
n_interior=length(interior_rho);
tip_z=s_fresnel(1);
center_z=d_ref;
n_ring_verts=n_theta*n_interior;
verts=zeros(n_ring_verts+2,3);
idx=1;
for it=1:n_theta
ct=cos(theta(it));st=sin(theta(it));
for jp=1:n_interior
verts(idx,:)=[interior_rho(jp)*ct,interior_rho(jp)*st,interior_z(jp)];
idx=idx+1;
end
end
tip_vertex_idx=n_ring_verts+1;
center_vertex_idx=n_ring_verts+2;
verts(tip_vertex_idx,:)=[0,0,tip_z];
verts(center_vertex_idx,:)=[0,0,center_z];
faces=zeros(n_theta*(n_interior-1)*2+n_theta*2,3);
fidx=1;
for it=1:n_theta
it_next=mod(it,n_theta)+1;
for jp=1:(n_interior-1)
a=(it-1)*n_interior+jp;
b=(it-1)*n_interior+jp+1;
cc=(it_next-1)*n_interior+jp+1;
d=(it_next-1)*n_interior+jp;
faces(fidx,:)=[a b cc];fidx=fidx+1;
faces(fidx,:)=[a cc d];fidx=fidx+1;
end
end
for it=1:n_theta
it_next=mod(it,n_theta)+1;
a=(it-1)*n_interior+1;
b=(it_next-1)*n_interior+1;
faces(fidx,:)=[a b tip_vertex_idx];fidx=fidx+1;
end
for it=1:n_theta
it_next=mod(it,n_theta)+1;
a=(it-1)*n_interior+n_interior;
b=(it_next-1)*n_interior+n_interior;
faces(fidx,:)=[a center_vertex_idx b];fidx=fidx+1;
end
fprintf('Mesh: %d vertices, %d faces\n',size(verts,1),size(faces,1));
stl_filename=sprintf('fresnel_lens_%.0fGHz_blank%.0fmm_active%.0fmm_f%.0fmm.stl',freq_GHz,diameter_mm,active_diameter_mm,f_work);
write_ascii_stl(stl_filename,verts,faces);
full_stl_path=fullfile(pwd,stl_filename);
if isfile(full_stl_path)
file_info=dir(full_stl_path);
fprintf('\nSTL successfully saved.\n');
fprintf(' File: %s\n',stl_filename);
fprintf(' Location: %s\n',full_stl_path);
fprintf(' Size: %.1f KB (%d triangles)\n',file_info.bytes/1024,size(faces,1));
else
error('STL write failed: %s',full_stl_path);
end
fprintf('Physical blank diameter: %.1f mm\n',2*blank_radius);
fprintf('Optically active diameter: %.1f mm\n',2*rho_max);
fprintf('Mounting rim width: %.2f mm\n',rim_width);
fprintf('Max thickness: %.2f mm, min thickness: %.2f mm\n',d_ref,d_ref-t_groove);

function val=valid_input(prompt_str,lo,hi)
while true
str_in=input(prompt_str,'s');
val=str2double(str_in);
if ~isnan(val)&&isscalar(val)&&val>=lo&&val<=hi
return;
else
fprintf('Invalid input. Enter a value between %.1f and %.1f.\n',lo,hi);
end
end
end

function write_ascii_stl(filename,verts,faces)
fid=fopen(filename,'w');
fprintf(fid,'solid fresnel_lens\n');
for i=1:size(faces,1)
v1=verts(faces(i,1),:);
v2=verts(faces(i,2),:);
v3=verts(faces(i,3),:);
nvec=cross(v2-v1,v3-v1);
nn=norm(nvec);
if nn>0
nvec=nvec/nn;
end
fprintf(fid,'facet normal %e %e %e\n',nvec(1),nvec(2),nvec(3));
fprintf(fid,' outer loop\n');
fprintf(fid,' vertex %e %e %e\n',v1(1),v1(2),v1(3));
fprintf(fid,' vertex %e %e %e\n',v2(1),v2(2),v2(3));
fprintf(fid,' vertex %e %e %e\n',v3(1),v3(2),v3(3));
fprintf(fid,' endloop\n');
fprintf(fid,'endfacet\n');
end
fprintf(fid,'endsolid fresnel_lens\n');
fclose(fid);
end