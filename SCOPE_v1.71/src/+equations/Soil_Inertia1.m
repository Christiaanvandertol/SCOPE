function [GAM] = Soil_Inertia1(SMC)

%soil inertia method by Murray and Verhoef (

%% parameters

theta_s = 0.435; %(saturated water content, m3/m3)
Sr = SMC/theta_s;

%fss = 0.58; %(sand fraction)
gamma_s = 0.96; %(soil texture dependent parameter)
dels = 1.33; %(shape parameter)


ke = exp(gamma_s*(1- power(Sr,(gamma_s - dels))));

phis  = 0.435; %(phis == theta_s)
lambda_d = -0.56*phis + 0.51;

QC = 0.60; %(quartz content)
lambda_qc = 7.7;  %(thermal conductivity of quartz, constant)

lambda_s = (lambda_qc^(QC))*lambda_d^(1-QC);
lambda_wtr = 0.57;   %(thermal conductivity of water, W/m.K, constant)

lambda_w = (lambda_s^(1-phis))*lambda_wtr^(phis);

lambdas = ke*(lambda_w - lambda_d) + lambda_d;

Hcs = 2*10^6;
Hcw = 4.2*10^6;

Hc = (Hcw * SMC)+ (1-theta_s)*Hcs;

GAM = sqrt(lambdas.*Hc);
