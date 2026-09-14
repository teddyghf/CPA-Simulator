'use strict';
// Public geometry/absorption, estimated LP01 mode properties. See physics-notes.html.
// beta2: ps^2/m; area: um^2; gamma: W^-1 m^-1. Pump controls are launched watts.
const FIBERS=Object.freeze({
 pre1:Object.freeze({name:'PM Yb 6/125',core:6,clad:125,na:.12,area:38.9135,beta2:.021652,absorption:250,pumpMode:'单模纤芯泵浦',maxPump:1,pumpKey:'ld1',lengthKey:'pre1Length',defaultLength:2}),
 pre2:Object.freeze({name:'PM Yb 14/125',core:14,clad:125,na:.07,area:144.0413,beta2:.018943,absorption:16.6,pumpMode:'包层泵浦',maxPump:9,pumpKey:'ld2',lengthKey:'pre2Length',defaultLength:2}),
 pre3:Object.freeze({name:'PM Yb 20/130',core:20,clad:130,na:.08,area:201.2571,beta2:.018477,absorption:10.2,pumpMode:'包层泵浦',maxPump:27,pumpKey:'ld3',lengthKey:'pre3Length',defaultLength:2}),
 coupler:Object.freeze({name:'aeroGAIN-ROD 3.1',core:85,clad:260,area:3300,beta2:.018973,absorption:17/.804,pumpMode:'反向包层泵浦',maxPump:300,pumpKey:'pump',defaultLength:.8})
});
const AMP_DEFAULTS=FIBERS;
const YB=Object.freeze({h:6.62607015e-34,c:299792458,lambda:1030e-9,pumpLambda:976e-9,n2:2.74e-20,lifetime:.001,pumpAbs:2.5e-24,pumpEm:2.5e-24,signalAbs:.07e-24,signalEm:.6e-24});
// Smooth reference cross-section curves, not a measured spectrum for a commercial batch.
function ybCrossSections(lambdaNm){const offset=lambdaNm-1030;return{a:YB.signalAbs*Math.exp(-offset/18),e:YB.signalEm*Math.exp(-4*Math.log(2)*(offset/45)**2)};}
const ASE_BINS=Array.from({length:80},(_,i)=>{const wavelength=1000+(i+.5)*1.25,cross=ybCrossSections(wavelength),frequency=YB.c/(wavelength*1e-9),bandwidth=YB.c*1.25e-9/(wavelength*1e-9)**2;return{wavelength,frequency,bandwidth,...cross};});
function fiberProperties(id,length){const f=FIBERS[id],coreArea=Math.PI*(f.core*1e-6/2)**2,overlap=1-Math.exp(-2*coreArea/(f.area*1e-12)),pumpOverlap=id==='pre1'?overlap:(f.core/f.clad)**2;return{...f,length:id==='coupler'?.8:(length??f.defaultLength),coreArea,overlap,pumpOverlap,density:f.absorption*Math.LN10/10/(pumpOverlap*YB.pumpAbs),gamma:2*Math.PI*YB.n2/(YB.lambda*f.area*1e-12)};}
function expIntegral(x){return Math.abs(x)<1e-6?1+x/2+x*x/6:Math.expm1(x)/x;}
// Uniform-inversion steady-state rate equation, including bidirectional ASE and pump depletion.
// One reservoir per gain fiber. Its integrated photon balance is solved, not an ASE percentage.
function solveAmplifier(id,pumpW,length,signalInW,aseIn=[]){
 const f=fiberProperties(id,length),L=f.length,N=f.density,gs=f.overlap*N,ep=YB.h*YB.c/YB.pumpLambda,es=YB.h*YB.c/YB.lambda,pump=Math.max(0,Math.min(f.maxPump,pumpW)),signal=Math.max(0,signalInW);
 const evaluate=inversion=>{
  const pumpOut=pump*Math.exp(-f.pumpOverlap*N*(YB.pumpAbs-(YB.pumpAbs+YB.pumpEm)*inversion)*L),logGain=gs*((YB.signalAbs+YB.signalEm)*inversion-YB.signalAbs)*L,gain=Math.exp(logGain),signalOut=signal*gain;
  let aseStimPhotons=0,aseW=0,backwardW=0;const bins=[];
  for(let i=0;i<ASE_BINS.length;i++){const b=ASE_BINS[i],x=gs*((b.a+b.e)*inversion-b.a)*L,ePhoton=YB.h*b.frequency,spontaneous=2*ePhoton*b.bandwidth*gs*b.e*inversion*L,created=spontaneous*expIntegral(x),input=aseIn[i]||0,forward=input*Math.exp(x)+created;
   aseStimPhotons+=(forward-input+created-2*spontaneous)/ePhoton;aseW+=forward;backwardW+=created;bins.push(forward);
  }
  const decayPhotons=inversion*N*f.coreArea*L/YB.lifetime,residual=(pump-pumpOut)/ep-(signalOut-signal)/es-aseStimPhotons-decayPhotons;
  return{inversion,pumpW:pump,pumpOut,absorbedPump:pump-pumpOut,signalIn:signal,signalOut,aseW,backwardW,aseBins:bins,gain,logGain,residual,decayPhotons,...f};
 };
 let lo=0,hi=.5;for(let i=0;i<60;i++){const mid=(lo+hi)/2;if(evaluate(mid).residual>0)lo=mid;else hi=mid;}
 return evaluate((lo+hi)/2);
}
function asePSDAt(bins,frequencyOffsetTHz){const wavelength=299792.458/(299792.458/1030+frequencyOffsetTHz),at=(wavelength-1000)/1.25-.5,i=Math.floor(at),fraction=at-i;if(i<0||i>=ASE_BINS.length-1)return 0;return ((bins[i]||0)/ASE_BINS[i].bandwidth*(1-fraction)+(bins[i+1]||0)/ASE_BINS[i+1].bandwidth*fraction)*1e12;}
function nonlinearIntegral(gamma,length,inputPeakW,gain){return gamma*length*inputPeakW*expIntegral(Math.log(Math.max(gain,1e-300)));}
