// FX estimation by SWIPE
//
// (c) 2015 Mark Huckvale University College London
//
if( 'undefined' === typeof window) importScripts("fft.js");

//
// function [p,t,s] = swipep(x,fs,plim,dt,dlog2p,dERBs,woverlap,sTHR)
// SWIPEP Pitch estimation using SWIPE.
//    P = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,STHR) estimates the pitch
//    of the vector signal X every DT seconds. The sampling frequency of
//    the signal is Fs (in Hertz). The spectrum is computed using a Hann
//    window with an overlap WOVERLAP between 0 and 1. The spectrum is
//    sampled uniformly in the ERB scale with a step size of DERBS ERBs. The
//    pitch is searched within the range [PMIN PMAX] (in Hertz) with samples
//    distributed every DLOG2P units on a base-2 logarithmic scale of Hertz.
//    The pitch is fine-tuned using parabolic interpolation with a resolution
//    of 1 cent. Pitch estimates with a strength lower than STHR are treated
//    as undefined.
//
//    [P,T,S] = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR)
//    returns the times T at which the pitch was estimated and the pitch
//    strength S of every pitch estimate.
//
//    P = SWIPEP(X,Fs) estimates the pitch using the default settings PMIN =
//    30 Hz, PMAX = 5000 Hz, DT = 0.001 s, DLOG2P = 1/48 (48 steps per
//    octave), DERBS = 0.1 ERBs, WOVERLAP = 0.5, and STHR = -Inf.
//
//    P = SWIPEP(X,Fs,...,[],...) uses the default setting for the parameter
//    replaced with the placeholder [].
//
//    REMARKS: (1) For better results, make DLOG2P and DERBS as small as
//    possible and WOVERLAP as large as possible. However, take into account
//    that the computational complexity of the algorithm is inversely
//    proportional to DLOG2P, DERBS and 1-WOVERLAP, and that the  default
//    values have been found empirically to produce good results. Consider
//    also that the computational complexity is directly proportional to the
//    number of octaves in the pitch search range, and therefore , it is
//    recommendable to restrict the search range to the expected range of
//    pitch, if any. (2) This code implements SWIPE, which uses only the
//    first and prime harmonics of the signal. To convert it into SWIPE,
//    which uses all the harmonics of the signal, replace the word
//    PRIMES with a colon (it is located almost at the end of the code).
//    However, this may not be recommendable since SWIPE is reported to
//    produce on average better results than SWIPE (Camacho and Harris,
//    2008).
//
//    REFERENCES: Camacho, A., Harris, J.G, (2008) A sawtooth waveform
//    inspired pitch estimator for speech and music, J. Acoust. Soc. Am.
//    124, 1638-1652.
//
//    MAINTENANCE HISTORY:
//    - Added line 153 to avoid division by zero in line 154 if loudness
//      equals zero (06/23/2010).

function log2(val) {
	return Math.log10(val)/Math.log10(2.0);
}
function hz2erbs(hz)
{
	return 6.44 * ( log2( 229 + hz ) - 7.84 );
}
function erbs2hz(erbs)
{
	return Math.pow(2,( erbs/6.44 + 7.84) ) - 229;
}

/* simple interpolation */
function interp1(x1,y1,x2)
{
	var y2=new Float32Array(x2.length);
	var n1=x1.length;
	var	i,j;
	var	x;
	var	m;

	for (i=0;i<x2.length;i++) {
		x=x2[i];
		if (x <= x1[0])
			y2[i] = y1[0];
		else if (x >= x1[n1-1])
			y2[i] = y1[n1-1];
		else {
			for (j=1;j<n1;j++) {
				if ((x1[j-1] <= x) && (x <= x1[j])) {
					m=(x1[j]-x)/(x1[j]-x1[j-1]);
					y2[i]=(m*y1[j-1]+(1-m)*y1[j]);
					break;
				}
			}
		}
	}
	return y2;
}

var plim=[40,600];
var dt=0.005;
var dlog2p=1.0/48.0;
var dERBs=0.1;
var woverlap=1.0/2.0;
var sTHR=-1.0E10;
var primes=[1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,447,53,59,61,67,71,73,79,83,89,97,
101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,
311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,
439,443,449,457,461,463,467,479,487,491,499];
var NPRIME=primes.length;
var maxnfERBs;

function allocvector(size)
{
	return new Float32Array(size);
}
function allocmatrix(nrow,ncol)
{
	var i;
//	console.log("allocmatrix: nrow="+nrow+" ncol="+ncol);
	var p=new Array(nrow);
	for (i=0;i<nrow;i++) p[i]=new Float32Array(ncol);
	return p;
}
function printvector(name,vec,len)
{
	var s=name+"="+vec[0];
	for (var i=1;i<len;i++) s += ","+vec[i];
	trace(s);
}

// -----------------------------------------------------------------------------------------
// Calculates the co-ordinates xe, ye of the extremum and returns
// the 2nd derivative d2 if the extremum is deemed reliable or 0.0 otherwise
// Negative return value indicates a valid maximum, positive return value a valid minimum.
// Zero return value indicates an unreliable extremum (this does not mean that d2 is zero).
// When d2 is exact zero, the function returns exact 0.0 in both xe and ye.
// http://www.ebyte.it/library/codesnippets/P3Interpolation.html
//
function P3Interpolate_Extremum(xl,xc,xu,yl,yc,yu)
{
	var	xe=xc,ye=yc;
	var	d1,d2;

	d2 = 2*((yu-yc)/(xu-xc)-(yl-yc)/(xl-xc))/(xu-xl);
	d1 = (yu-yc)/(xu-xc) - 0.5*d2*(xu-xc);
	if (d2) {
		xe = (xc - d1/d2);
		ye = (yc + 0.5*d1*(xe-xc));
	}
	if ((xe<xl)||(xe>xu)) {
		// fail
		xe=xc;
		ye=yc;
	}

	return { xe:xe, ye:ye };
}

var k=[],q=[],c2pq=[];

function pitchStrengthOneCandidate(f,foffset,NL,nf,pc)
{
	var	n=Math.floor(f[foffset+nf-1]/pc-0.75);		// number of harmonics
	var	i,j;
	var	sum,norm;

	if (k.length==0) k=allocvector(maxnfERBs);
	if (q.length==0) q=allocvector(maxnfERBs);
	if (c2pq.length==0) c2pq=allocvector(maxnfERBs);

	for (i=0;i<nf;i++) {
		k[i]=0;
		q[i]=f[foffset+i]/pc;
		c2pq[i]=Math.cos(2*Math.PI*q[i]);
	}

	for (i=0;(i<n) && (i<NPRIME);i++) {
		// binary search for start sweep position
		var ii=0;
		var kk=nf-1;
		var jj;
		do {
			jj=Math.floor((ii+kk)/2);
			if (q[jj]<primes[i])
				ii=jj+1;
			else if (q[jj]>primes[i])
				kk=jj-1;
			else
				break;
		} while (ii < kk);
		while ((jj>0)&&(q[jj]>(primes[i]-0.75))) jj--;
		// sweep for ratios close to prime target
		for (j=jj;j<nf;j++) {
			var a=(q[j]-primes[i]);
			if (a < -0.75)
				/* skip */;
			else if (a < -0.25)
				k[j] += (c2pq[j]/2);
			else if (a < 0.25)
				k[j] += c2pq[j];
			else if (a < 0.75)
				k[j] += (c2pq[j]/2);
			else
				break;
		}
	}
	sum=0;
	for (i=0;i<nf;i++) {
		k[i] = (k[i]*Math.sqrt(1/f[foffset+i]));
		if (k[i] > 0) sum+=k[i]*k[i];
	}
	norm=Math.sqrt(sum);

	sum=0;
	for (i=0;i<nf;i++) sum += k[i] * NL[i] / norm;

	return(sum);
}



function SWIPEP(x,fs,fxlen)
{
	var	xlen=x.length;
	var	fx=allocvector(fxlen);
	var vs=allocvector(fxlen);
	var	i,j,k,l;
	var	t;
	var	nt;
	var	log2pc,pc;
	var	nlog2pc,npc;
	var	S;
	var	logWs=new Array(2);
	var	ws;
	var	p0;
	var	nws;
	var	d;
	var	nd;
	var	dn;
	var	xzp;
	var	w;
	var	omega;
	var	o;
	var	X;
	var	nX,nfft;
	var	f,ti;
	var	fERBs;
	var	nfERBs;
	var	npad;
	var	jcand,kcand;
	var	njcand,nkcand;
	var	L;
	var	Si;
	var	N;
	var	total;
	var	kint;
	var	nkint;
	var	norm;
	var	NL;
	var	lambda,mu;
	var	smax;
	var	jmax;
	var	fint,sint;
	var	sold,snew;

	// calculate window times
	nt=1+Math.floor((xlen/fs)/dt);
	t=allocvector(nt);
	for (i=0;i<nt;i++) t[i]=(i*dt);
//	printvector("t",t,nt);

	// Define pitch candidates
	nlog2pc=1+Math.floor((log2(plim[1])-log2(plim[0]))/dlog2p);
	log2pc=allocvector(nlog2pc);
	for (i=0;i<nlog2pc;i++)
		log2pc[i] = (log2(plim[0])+i*dlog2p);
//	printvector("log2pc",log2pc,nlog2pc);

	npc=nlog2pc;
	pc=allocvector(npc);
	for (i=0;i<npc;i++)
		pc[i] = (Math.pow(2.0,log2pc[i]));
//	printvector("pc",pc,npc);

	// pitch strength matrix
	S=allocmatrix(nt,npc);

	// Determine P2-WSs
	logWs[0] = Math.floor(0.5+( log2( 8*fs / plim[0] ) ));
	logWs[1] = Math.floor(0.5+( log2( 8*fs / plim[1] ) ));
	nws = 1+Math.floor(logWs[0]-logWs[1]);
//	printf("logws[0]=%d logws[1]=%d nws=%d\n",logWs[0],logWs[1],nws);
	ws=new Int32Array(nws);
	for (i=0;i<nws;i++)
		ws[i] = Math.pow(2.0, logWs[0]-i);
//	printvector("ws",ws,nws);
	p0 = allocvector(nws);
	for (i=0;i<nws;i++)
		p0[i] = (8.0 * fs / ws[i]);	// Optimal pitches for P2-WSs
//	printvector("p0",p0,nws);

	// Determine window sizes used by each pitch candidate
	nd=nlog2pc;
	d=allocvector(nd);
	for (i=0;i<nd;i++)
		d[i] = (1 + log2pc[i] - log2( 8.0*fs/ws[0] ));
//	printvector("d",d,nd);

	// Create ERB-scale uniformly-spaced frequencies (in Hertz)
	nfERBs = 1+Math.floor((hz2erbs(fs/2)-hz2erbs(pc[0]/4))/dERBs);
	maxnfERBs=nfERBs;
	fERBs = allocvector(nfERBs);
	for (i=0;i<nfERBs;i++)
		fERBs[i] = erbs2hz(hz2erbs(pc[0]/4)+i*dERBs);
//	printvector("fERBs",fERBs,nfERBs);

	var fft = new RealFFT();

	// main loop over pitch candidates
	for (i=0;i<nws;i++) {

//		printf("i=%d nws=%d ws[i]=%d p0[i]=%g npc=%d\n",i,nws,ws[i],p0[i],npc);

		//calculate hop size
		dn = Math.floor(0.5+8*(1-woverlap)*fs/p0[i]);
		if (dn < 1) dn=1;
//		printf("dn=%d\n",dn);

		// Zero pad signal
		npad = ws[i]/2;
		xzp=allocvector(xlen+3*npad);
		for (j=0;j<xlen;j++) xzp[npad+j]=x[j];
//		printf("npad=%d xlen=%d xzplen=%d\n",npad,xlen,xlen+3*npad);

	    // Compute Hann window
	    w=allocvector(ws[i]);
	    omega=2*Math.PI/(ws[i]-1);
	    for (j=0;j<ws[i];j++) w[j]=0.5-0.5*Math.cos(j*omega);

		// calculate overlap
		o = ws[i]-dn;
		if (o < 0) o=0;
//		printf("ws[i]=%d dn=%d o=%d\n",ws[i],dn,o);


		// compute spectra
		nX=Math.floor((xlen+2*npad)/(ws[i]-o));
		nfft=2;
		while (nfft < ws[i]) nfft *= 2;
//		printf("nX=%d nfft=%d\n",nX,nfft);
		X=allocmatrix(nX,nfft);
		for (j=0;j<nX;j++) {
			for (k=0;k<ws[i];k++)
				X[j][k] = xzp[j*(ws[i]-o)+k] * w[k];

			X[j]=fft.RealFFTPower(X[j],ws[i]);
//			trace("ws[i]="+ws[i]+" nfft="+nfft+" X[j].length="+X[j].length);
//			X[j][1]=0;

//			for (k=0;k<nfft/2;k++) X[j][k] = sqrt(X[j][2*k]*X[j][2*k]+X[j][2*k+1]*X[j][2*k+1]);
		//	printvector("X[j]",X[j],nfft/2);
		}

		// FFT frequencies
		f=allocvector(nfft/2);
		for (k=0;k<nfft/2;k++) f[k] = (k*fs/nfft);
		//printvector("f",f,nfft/2);
		ti=allocvector(nX);
		for (k=0;k<nX;k++) ti[k] = (k*(ws[i]-o)/fs);
		//printvector("ti",ti,nX);

		// Select candidates that use this window size
		jcand=new Int32Array(nlog2pc);
		kcand=new Int32Array(nlog2pc);
		njcand=0;
		nkcand=0;
		if (nws==1) {
			for (k=0;k<nlog2pc;k++) { jcand[njcand++]=k; }
			nkcand=0;
		}
		else if ((i+1)==nws) {
			for (k=0;k<nlog2pc;k++) if ((d[k]-i-1)>-1) { jcand[njcand++]=k; }
			for (k=0;k<njcand;k++) if ((d[jcand[k]]-i-1)<0) { kcand[nkcand++]=k; }
		}
		else if (i==0) {
			for (k=0;k<nlog2pc;k++) if ((d[k]-1)<1) { jcand[njcand++]=k; }
			for (k=0;k<njcand;k++) if ((d[jcand[k]]-1)>0) { kcand[nkcand++]=k; }
		}
		else {
			for (k=0;k<nlog2pc;k++) if (Math.abs(d[k]-i-1)<1) { jcand[njcand++]=k; }
			for (k=0;k<njcand;k++) kcand[nkcand++]=k;
		}
		// printvector("jcand",jcand,njcand);
		// printvector("kcand",kcand,nkcand);

		// Compute loudness at ERBs uniformly-spaced frequencies
		for (j=0;j<nfERBs;j++)
			if (fERBs[j] > pc[jcand[0]]/4) break;
		//printf("npc=%d nfERBs=%d j=%d\n",npc,nfERBs,j);
		//fflush(stdout);
		for (k=0;k<nfERBs-j;k++) fERBs[k]=fERBs[j+k];
		nfERBs=nfERBs-j;
		//printvector("fERBs",fERBs,nfERBs);
		//fflush(stdout);

		L=allocmatrix(nX,nfERBs);
		for (j=0;j<nX;j++) {
			L[j]=interp1(f,X[j],fERBs);
			for (k=0;k<nfERBs;k++)
				L[j][k] = Math.sqrt((L[j][k]>0) ? L[j][k] : 0);
		//	printvector("L[j]",L[j],nfERBs);
		}

		// integration matrix
		nkint=njcand+1;
		kint=new Int32Array(nkint);
		kint[0]=0;
		for (j=0;j<nkint-1;j++) {
			for (k=0;k<nfERBs-kint[j];k++) {
				if (fERBs[kint[j]+k] > pc[jcand[j]]/4) break;
			}
			kint[j+1]=kint[j]+k;
		}
		for (j=0;j<nkint-1;j++) kint[j]=kint[j+1]+1;
		nkint--;
//		printivector("kint",kint,nkint);

		// Create loudness normalization matrix
		N=allocmatrix(nX,nfERBs);
		for (j=0;j<nX;j++) {
			total=0;
			for (k=nfERBs-1;k>=0;k--) {
				total += L[j][k]*L[j][k];
				N[j][k] = Math.sqrt(total);
			}
		//	printvector("N[j]",N[j],nfERBs);
		}

		// pitch strength matrix
		Si=allocmatrix(nX,njcand);

		// normalised loudness
		NL=allocmatrix(nX,nfERBs);

		for (j=0;j<njcand;j++) {
			for (k=0;k<nX;k++) {
				norm=N[k][kint[j]];
//				printf("j=%d k=%d kint[j]=%d norm=%g\n",j,k,kint[j],norm);
				for (l=kint[j]-1;l<nfERBs;l++) {
					NL[k][l-kint[j]+1]=L[k][l]/norm;
				}
				// printvector("NL[k]",NL[k],nfERBs-kint[j]+1);
				Si[k][j]=pitchStrengthOneCandidate(fERBs,kint[j]-1,NL[k],nfERBs-kint[j]+1,pc[jcand[j]]);
			}
		}
//		for (k=0;k<nX;k++) printvector("Si[k]",Si[k],njcand);

		lambda = allocvector(njcand);
		for (j=0;j<njcand;j++) lambda[j] = d[jcand[kcand[j]]]-i-1;
//		printvector("lambda",lambda,njcand);
		mu = allocvector(njcand);
		for (j=0;j<njcand;j++) mu[j] = 1-Math.abs(lambda[j]);
//		printvector("mu",mu,njcand);
//		fflush(stdout);

		// add into accumulator
		sold=allocvector(nX);
		snew=allocvector(nt);
		for (j=0;j<njcand;j++) {
			// interpolate Si to required times
			for (k=0;k<nX;k++) sold[k]=Si[k][j];
			snew=interp1(ti,sold,t);
			// accumulate
			for (k=0;k<nt;k++) {
				S[k][jcand[j]] += mu[j]*snew[k];
			}
		}
//		for (k=0;k<nt;k++) printvector("S[k]",S[k],npc);


	}

//	printf("Done.\n");  fflush(stdout);

	// find candidate with maximum value in each frame
	for (i=0;i<nt && i<fxlen;i++) {
		smax=sTHR;
		jmax=0;
		for (j=0;j<npc;j++) {
			if (S[i][j] > smax) {
				smax = S[i][j];
				jmax=j;
			}
		}
//		printf("i=%d smax=%g jmax=%d npc=%d nt=%d fxlen=%d\n",i,smax,jmax,npc,nt,fxlen);
//		fflush(stdout);
		if (smax <= sTHR) {
			fx[i]=0;
			vs[i]=smax;
		}
		else if ((jmax==0)||(jmax==npc-1)) {
			fx[i]=pc[jmax];
			vs[i]=smax;
		}
		else {
			// parabolic interpolation
			var val=P3Interpolate_Extremum(pc[jmax-1],pc[jmax],pc[jmax+1],S[i][jmax-1],S[i][jmax],S[i][jmax+1]);
//			printf("%g,%g,%g and %g,%g,%g goto %g,%g\n",pc[jmax-1],pc[jmax],pc[jmax+1],S[i][jmax-1],S[i][jmax],S[i][jmax+1],fint,sint);
			fx[i]=val.xe;
			vs[i]=val.ye;
		}
	}


	return { fx:fx, vs:vs, fxlen:((nt<fxlen)?(nt):(fxlen)) };
}


// do signal calculations
self.onmessage=function(event)
{
	var signal=event.data.signal;		// input signal
	var siglen=signal.length;		// input signal
	var srate=event.data.srate;
 	var fxlen=Math.round((siglen/srate)*200);

	// calculate FX
	var res = SWIPEP(signal,srate,fxlen);	
	
	// send everything back
	postMessage(res);

}

