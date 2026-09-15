'use strict';
const clamp=(x,lo,hi)=>Math.max(lo,Math.min(hi,x));
// Teaching AOM calibration: exact 0 V -> 0%, 1 V -> 90%; flat endpoints.
function aomEfficiency(voltage){const v=clamp(voltage,0,1);return .9*(6*v**5-15*v**4+10*v**3);}
function pidAdvance(s,available,dt){
 const measured=available*s.aomT/20,e=s.target/20-measured;
 const slope=-(measured-s.pidPrevious)/dt;
 s.pidDerivative+=(slope-s.pidDerivative)*(1-Math.exp(-dt/.06));s.pidPrevious=measured;
 let command=s.manualVoltage;
 if(s.powerLock){const p=s.kp*e,d=s.kd*s.pidDerivative,raw=p+s.integral+d;
  if((raw>0&&raw<1)||(raw>=1&&e<0)||(raw<=0&&e>0))s.integral=clamp(s.integral+s.ki*e*dt,-20,20);
  command=clamp(p+s.integral+d,0,1);s.pidTerms={p,i:s.integral,d,command};
 }else s.pidTerms={p:0,i:0,d:0,command};
 // Fixed 150 ms RF/control latency and 120 ms actuator inertia, slowed for teaching.
 s.voltageQueue??=Array(15).fill(s.aomV||0);s.voltageQueue.push(command);
 const delayed=s.voltageQueue.shift();s.aomV=clamp(s.aomV+(delayed-s.aomV)*(1-Math.exp(-dt/.12)),0,1);
 s.aomT=aomEfficiency(s.aomV);return available*s.aomT;
}
const MOTOR={min:0,max:50,optimum:12,mmPerChirp:2,baseDistance:80};
function compressorSample(position,tau0){const chirp=(position-MOTOR.optimum)/MOTOR.mmPerChirp,tau=tau0*Math.sqrt(1+chirp**2);return{position,distance:MOTOR.baseDistance+position,chirp,tau,signal:tau0/tau};}
function signalLoss(signal,targetSignal){return (1/Math.max(signal,1e-9)-1/targetSignal)**2;}
// The controller knows only sampled PD5 signals and a target calibration, never the plant optimum.
function pd5GradientStep(position,readSignal,targetSignal,learningRate,motorStep){
 const left=clamp(position-.05,MOTOR.min,MOTOR.max),right=clamp(position+.05,MOTOR.min,MOTOR.max);
 const signalLeft=readSignal(left),signalRight=readSignal(right);
 const gradient=(signalLoss(signalRight,targetSignal)-signalLoss(signalLeft,targetSignal))/Math.max(right-left,1e-9);
 const next=clamp(position-learningRate*motorStep*gradient,MOTOR.min,MOTOR.max);
 return{position:next,gradient,left,right,signalLeft,signalRight,move:next-position};
}
function compressionReachable(tau0,target){return target>=tau0&&target<=Math.max(compressorSample(MOTOR.min,tau0).tau,compressorSample(MOTOR.max,tau0).tau);}
function gaussianOutput(tau0,tau){
 const ln2=Math.log(2),chirp=Math.sqrt(Math.max(0,(tau/tau0)**2-1)),gdd=chirp*tau0*tau0/(4*ln2),bandwidth=2*ln2/Math.PI*1000/tau0;
 const temporal=[],spectral=[];
 for(let i=0;i<=240;i++){const z=-1.5+3*i/240,t=z*tau,nu=z*bandwidth,intensity=Math.exp(-4*ln2*z*z);
  temporal.push({x:t/1000,intensity,phase:-2*ln2*chirp*z*z});
  spectral.push({x:nu,intensity,phase:.5*gdd*(2*Math.PI*nu/1000)**2});
 }
 return{tau0,tau,chirp,gdd,bandwidth,temporal,spectral};
}