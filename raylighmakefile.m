clc
clear all
L=[64, 256];

rhostar = [3,1000];
mustar  = [1,100];

At = [0.5,    0.998];
Re = [3000,   3000];
Ca = [0.26,   0.44]; 
Pe = [1000,   1000];
to = [100,  8000];

W = [5,       5];

rho = [3,     1];
rhol= rho ./ rhostar;

g  = L./(to.^2.*At);
Uo = sqrt(L.*g);
muh= rho.*Uo.*L./Re;
mul= muh./mustar;
nuh=muh./rho;
nul=mul./(rho./rhostar);

sigma = muh.*Uo./Ca;
M = Uo.*L ./ Pe;

names = ['RTI_LDR.xml', 'RTI_HDR.xml'];
time = [3, 2];

% for (i in [1,2)){
%     print(names[i])
%     sink(names[i]) ?>
% <?xml version="1.0"?>
% <!--Model:	d2q9_pf_velocity 
%     Created: 	17-01-2020 
%     By:		T.Mitchell -->
% <CLBConfig version="2.0" output="output./" permissive="true">
% 	<Geometry nx="<?%f L[i] ?>" ny="<?%f 4.*L[i] ?>">
% 		<MRT>
% 			<Box./>
% 		<./MRT>
% 		<SpikeTrack>
% 			<Box dx="<?%f 0.5.*L[i] ?>" nx="1" dy="10" fy="<?%f 4.*L[i] - 10 ?>"./>
% 		<./SpikeTrack>
% 		<BubbleTrack>
% 			<Box nx="1" dy="10" fy="<?%f 4.*L[i] - 10 ?>"./>
% 		<./BubbleTrack>
% 
% 		<Wall mask="ALL">
% 			<Box dy="-1"./>
% 		<./Wall>
% 		<Wall mask="ALL">
% 			<Box ny="1"./>
% 		<./Wall>
% 	<./Geometry>
% 	<Model>
% 		<Param name="MidPoint" value="<?%f 2.*L[i] ?>"./>
% 		<Param name="Perturbation" value="0.1"./>
% 		<Param name="Period" value="<?%f L[i] ?>"./>
% 		<Param name="sigma" value="<?%f sigma[i] ?>"./>
% 		<Param name="M" value="<?%f M[i] ?>"./>
% 		<Param name="W" value="<?%f W[i] ?>"./>
% 		<Param name="GravitationY" value="<?%.8f -1.*g[i] ?>"./>
% 		<Param name="Viscosity_h" value="<?%.8f nuh[i] ?>"./>
% 		<Param name="Viscosity_l" value="<?%.8f nul[i] ?>"./>
% 		<Param name="PhaseField_init" value="1.0"./>
% 		<Param name="Density_l" value="<?%f rhol[i] ?>"./>
% 		<Param name="Density_h" value="<?%f rho[i]  ?>"./>
% 	<./Model>
% 	<VTK./>
% 	<Solve Iterations="<?%f time[i] .* to[i] ?>" output="output./">
% 		<VTK Iterations="<?%f  0.5.*to[i] ?>"./>
% 		<Log Iterations="<?%f 0.01.*to[i] ?>" ./>
% 	<./Solve>
% <./CLBConfig>
% <?R
%     sink()
% }
% ?>
