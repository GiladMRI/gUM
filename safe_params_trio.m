function X=safe_params_trio()
X=struct('name','trio');
X.look_ahead=1.02; % ???
% Online tolerance 1 0.22, 2 0.4872, Correction factor XYZ 0.98
% stimulus correction factor 1, SW correction factor 1
% T limit 0.0451128, 3
% delta T 0.01, Limits 20,1330, Dbdt threshold 20
X.x.tau1=0.131;
X.x.tau2=12;
X.x.tau3=0.5;
X.x.a1=0.3275;
X.x.a2=0.24;
X.x.a3=0.4325;
X.x.stim_limit=45.7371;
X.x.stim_thresh=36.5897;

X.y.tau1=0.1775;
X.y.tau2=12;
X.y.tau3=0.703;
X.y.a1=0.3265;
X.y.a2=0.24;
X.y.a3=0.4335;
X.y.stim_limit=27.6493;
X.y.stim_thresh=22.1194;

X.z.tau1=0.15;
X.z.tau2=12;
X.z.tau3=0.549;
X.z.a1=0.282;
X.z.a2=0.24;
X.z.a3=0.478;
X.z.stim_limit=31.9437;
X.z.stim_thresh=25.555;