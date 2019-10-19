/* filter -- General Purpose IIR Digital Filter Routines */

/* M.A.Huckvale -- November 1991 */

/* Code borrowed largely and loosely from:
		Cappellini, Constantinides, Emillani,
		"Digital Filters and their Applications"
		Academic Press, 1976
*/

function trace(message){
	if (typeof console == "object"){
		console.log(message);
	}
}

Complex=function() {
	this.re=0;
	this.im=0;
};

FilterSection=function() {
	this.ncoeff=0;
	this.acoeff=new Array(5);
	this.bcoeff=new Array(5);
	this.memory=[0,0,0,0,0,0];
};

Filter=function() {
	var self=this;
	var me=this.constructor.prototype;

	/* maximum sizes */
	this.MAX_ORDER=20;

	/* filter types */
	this.LOW_PASS=1;
	this.HIGH_PASS=2;
	this.BAND_PASS=3;
	this.BAND_STOP=4;

	/* structure of the filter */
	this.nsection=0;
	this.section=0;		/* array of FilterSection */

	/* filter design */
	this.a=new Array(this.MAX_ORDER);
	this.b=new Array(this.MAX_ORDER);

	/* do filter transformation */
	me.transform=function(nsection,fp,f1,f2,ftype)
	{
		var	al,a0,a1,a2,b0,b1,b2;
		var ncoeff;
		var cappa;
		var an=new Array(3);
		var fina=new Array(3);
		var finb=new Array(3);
		var	i,j,k;

		switch (ftype) {
		case this.LOW_PASS:
			al = Math.sin(Math.PI*(fp-f1))/Math.sin(Math.PI*(fp+f1));
			a0 = -al;
			a1 = 1.0;
			a2 = 0.0;
			b0 = a1;
			b1 = a0;
			b2 = 0.0;
			ncoeff=3;
			break;
		case this.HIGH_PASS:
			al = -Math.cos(Math.PI*(fp+f1))/Math.cos(Math.PI*(fp-f1));
			a0 = -al;
			a1 = -1.0;
			a2 = 0.0;
			b0 = -a1;
			b1 = -a0;
			b2 = 0.0;
			ncoeff = 3;
			break;
		case this.BAND_PASS:
			al = Math.cos(Math.PI*(f2+f1))/Math.cos(Math.PI*(f2-f1));
			cappa = Math.tan(Math.PI*fp)/Math.tan(Math.PI*(f2-f1));
			a0 = -(cappa-1.0)/(cappa+1.0);
			a1 = 2.0*al*cappa/(cappa+1.0);
			a2 = -1.0;
			b0 = -a2;
			b1 = -a1;
			b2 = -a0;
			ncoeff = 5;
			break;
		case this.BAND_STOP:
			al = Math.cos(Math.PI*(f2+f1))/Math.cos(Math.PI*(f2-f1));
			cappa = Math.tan(Math.PI*fp)*Math.tan(Math.PI*(f2-f1));
			a0 = (1.0-cappa)/(cappa+1.0);
			a1 = -2.0*al/(cappa+1.0);
			a2 = 1.0;
			b0 = a2;
			b1 = a1;
			b2 = a0;
			ncoeff = 5;
			break;
		default:
			throw "filter: unknown filter type: " + ftype;
			return(3);
		}

		an[0]=new Array(5);
		an[0][0] = b0*b0;
		an[0][1] = 2.0*b0*b1;
		an[0][2] = b1*b1 + 2.0*b0*b2;
		an[0][3] = 2.0*b1*b2;
		an[0][4] = b2*b2;

		an[1]=new Array(5);
		an[1][0] = a0*b0;
		an[1][1] = a0*b1+a1*b0;
		an[1][2] = a0*b2 + a1*b1 + a2*b0;
		an[1][3] = a1*b2 + a2*b1;
		an[1][4] = a2*b2;

		an[2]=new Array(5);
		an[2][0] = a0*a0;
		an[2][1] = 2.0 * a0 * a1;
		an[2][2] = a1*a1 + 2.0*a0*a2;
		an[2][3] = 2.0*a1*a2;
		an[2][4] = a2*a2;

		for (i=0;i<nsection;i++) {
			for (j=0;j<3;j++) {
				fina[j] = this.a[i][j];
				finb[j] = this.b[i][j];
			}
			for (j=0;j<5;j++) {
				this.a[i][j] = 0.0;
				this.b[i][j] = 0.0;
				for (k=0;k<3;k++) {
					this.a[i][j] += an[k][j] * fina[k];
					this.b[i][j] += an[k][j] * finb[k];
				}
			}
		}
		for (i=0;i<nsection;i++) {
			for (j=1;j<ncoeff;j++) {
				this.a[i][j] /= this.b[i][0];
				this.b[i][j] /= this.b[i][0];
			}
			this.a[i][0] /= this.b[i][0];
			this.b[i][0] = 1.0;
		}

		return(ncoeff);
	}

	/* filter design program */
	me.design=function(ftype,forder,flow,fhigh,fsamp)
	{
		var polre=new Array(this.MAX_ORDER);
		var polim=new Array(this.MAX_ORDER);
		var	ftrans;		/* freq, prior to transform */
		var	fca;
		var	i,j;
		var	ncoeff;

		switch (ftype) {
		case this.LOW_PASS:
			ftrans = flow/fsamp;
			break;
		case this.HIGH_PASS:
			ftrans = (fsamp-flow)/fsamp;
			break;
		case this.BAND_PASS:
			ftrans = (fhigh-flow)/fsamp;
			break;
		case this.BAND_STOP:
			ftrans = (fhigh-flow)/fsamp;
			break;
		default:
			throw "filter: unknown filter type: " +  ftype;
			return;
		}

		if (forder > this.MAX_ORDER) {
			throw "filter: order of filter too large"
			return;
		}
		if (forder & 1) {
			throw "filter: filter order must be even"
			return;
		}

		/* get filter corner frequency */
		fca = Math.tan(Math.PI*ftrans);

		/* determine pole and zero locations of prototype */
		this.nsection = forder/2;
		for (i=0;i<this.nsection;i++) {
			var tetap;
			tetap = Math.PI * (2.0 * i + 1) / (2.0 * forder);
			polre[i] = -fca * Math.cos(tetap);
			polim[i] = fca * Math.sin(tetap);
		}
		for (i=0;i<this.nsection;i++) {
			var p1,p2,p3;
			p1 = polim[i] * polim[i];
			p2 = (1.0 - polre[i]) * (1.0 - polre[i]);
			p3 = polre[i] * polre[i];
			polre[i] = (1.0 - p3 - p1)/(p2+p1);
			polim[i] = (2.0 * polim[i])/(p2+p1);
		}

		/* calculate coefficients from pole/zero spec */
		for (i=0;i<this.nsection;i++) {
			this.b[i]=new Array(5);
			this.b[i][0] = 1.0;
			this.b[i][1] = -2.0 * polre[i];
			this.b[i][2] = polre[i] * polre[i] + polim[i] * polim[i];
			this.b[i][3] = 0;
			this.b[i][4] = 0;
		}
		for (i=0;i<this.nsection;i++) {
			var tot;
			tot = 4.0/(this.b[i][0] + this.b[i][1] + this.b[i][2]);
			this.a[i]=new Array(5);
			this.a[i][0] = 1.0/tot;
			this.a[i][1] = 2.0/tot;
			this.a[i][2] = 1.0/tot;
			this.a[i][3] = 0;
			this.a[i][4] = 0;
		}

		if (ftype != this.LOW_PASS)
			/* transform low-pass, using a[],b[] global coeff */
			ncoeff = this.transform(this.nsection,ftrans,flow/fsamp,fhigh/fsamp,ftype);
		else
			ncoeff = 3;

		/* build FILTER structure */
		this.section=new Array(this.nsection);
		for (i=0;i<this.nsection;i++) {
			this.section[i]=new FilterSection();
			this.section[i].ncoeff = ncoeff;
			for (j=0;j<ncoeff;j++) {
				this.section[i].acoeff[j] = this.a[i][j];
				this.section[i].bcoeff[j] = this.b[i][j];
			}
		}
	}

	/* build a simple resonator */
	me.resonator=function(cf,bw,fsamp)
	{
		var	wf,wb,r,gain;

		/* build FILTER structure */
		this.nsection=1;
		this.section=new Array(1);
		this.section[0]=new FilterSection();
		this.section[0].ncoeff = 3;

		wf = 2*Math.PI*cf/fsamp;
		wb = 2*Math.PI*bw/fsamp;
		r = 1.0 - wb/2;

		this.section[0].acoeff[0] = 1.0;
		this.section[0].acoeff[1] = 0.0;
		this.section[0].acoeff[2] = 0.0; // -1.0;
		this.section[0].bcoeff[0] = 1.0;
		this.section[0].bcoeff[1] = -2.0 * r * Math.cos(wf);
		this.section[0].bcoeff[2] = r * r;

		/* set gain at DC to unity */
		gain = this.response(0,fsamp);
		this.section[0].acoeff[0] /= gain;
		this.section[0].acoeff[2] /= gain;

	}
	
	/* build an amplifier */
	me.amplifier=function(gain)
	{
		this.nsection=1;
		this.section=new Array(1);
		this.section[0]=new FilterSection();
		this.section[0].ncoeff = 1;
		this.section[0].acoeff[0] = gain;
		this.section[0].bcoeff[0] = 1.0;
	}

	/* build a vocal tract */
	me.vocaltract=function(f1,f2,f3,fsamp)
	{
		var i;
		var	wf,wb,r,gain;

		/* build FILTER structure */
		this.nsection=3;
		this.section=new Array(3);
		for (i=0;i<3;i++) {
			this.section[i]=new FilterSection();
			this.section[i].ncoeff = 3;
			this.section[i].acoeff[0] = 1.0;
			this.section[i].acoeff[1] = 0.0;
			this.section[i].acoeff[2] = 0.0; // -1.0;
			this.section[i].bcoeff[0] = 1.0;
		}
	
		wf = 2*Math.PI*f1/fsamp;
		wb = 2*Math.PI*60/fsamp;
		r = 1.0 - wb/2;
		this.section[0].bcoeff[1] = -2.0 * r * Math.cos(wf);
		this.section[0].bcoeff[2] = r * r;

		wf = 2*Math.PI*f2/fsamp;
		wb = 2*Math.PI*90/fsamp;
		r = 1.0 - wb/2;
		this.section[1].bcoeff[1] = -2.0 * r * Math.cos(wf);
		this.section[1].bcoeff[2] = r * r;

		wf = 2*Math.PI*f3/fsamp;
		wb = 2*Math.PI*120/fsamp;
		r = 1.0 - wb/2;
		this.section[2].bcoeff[1] = -2.0 * r * Math.cos(wf);
		this.section[2].bcoeff[2] = r * r;

		/* set gain at DC to unity */
		gain = this.response(0,fsamp);
		this.section[0].acoeff[0] /= gain;
		
	}
	
	/* zero filter memory */
	me.clear=function()
	{
		var i,j;

		for (i=0;i<this.nsection;i++) {
			for (j=0;j<=this.section[i].ncoeff;j++) {
				this.section[i].memory[j]=0;
			}
		}
	}

	/* filter a sample */
	me.sample=function(y)
	{
		var i,j;
		var mm,mmm;

		for (i=0;i<this.nsection;i++) {
			mm = this.section[i].ncoeff-1;
			this.section[i].memory[mm] = y;
			y = 0.0;
			for (j=0;j<this.section[i].ncoeff-1;j++) {
				mmm = mm - j;
				this.section[i].memory[mm] -= this.section[i].bcoeff[mmm]*this.section[i].memory[j];
				y += this.section[i].acoeff[mmm] * this.section[i].memory[j];
				this.section[i].memory[j] = this.section[i].memory[j+1];
			}
			y += this.section[i].acoeff[0] * this.section[i].memory[mm];
		}
		return(y);
	}

	/* filter a whole signal */
	me.process=function(signal)
	{
		var	i;
		var osignal=new Float32Array(new ArrayBuffer(4*signal.length));
		this.clear();
		for (i=0;i<signal.length;i++)
			osignal[i] = this.sample(signal[i]);
		return(osignal);		
	}

	/* complex multiply */
	me.CMLT=function(a,b)
	{
		var t = new Complex();
		t.re = a.re*b.re - a.im*b.im;
		t.im = a.re*b.im + a.im*b.re;
		return(t);
	}

	/* calculate response of filter at given frequency */
	me.response=function(f,fsamp)
	{
		var	i,j;
		var	omega;
		var	gain;
		var wj=new Array(5);
		var h=new Complex();
		var an=new Complex();
		var ad=new Complex();
		var	t=new Complex();
		var	r1,r2;

		omega = 2.0 * Math.PI * f/fsamp;
		for (i=0;i<5;i++) wj[i]=new Complex();
		wj[0].re = 1.0;
		wj[0].im = 0.0;
		wj[1].re = Math.cos(omega);
		wj[1].im = Math.sin(omega);
		wj[2]=this.CMLT(wj[1],wj[1]);
		wj[3]=this.CMLT(wj[2],wj[1]);
		wj[4]=this.CMLT(wj[2],wj[2]);

		h.re = 1.0;
		h.im = 0.0;
		for (i=0;i<this.nsection;i++) {
			an.re = 0;
			an.im = 0;
			ad.re = 0;
			ad.im = 0;
			for (j=0;j<this.section[i].ncoeff;j++) {
				an.re += wj[j].re * this.section[i].acoeff[j];
				an.im += wj[j].im * this.section[i].acoeff[j];
				ad.re += wj[j].re * this.section[i].bcoeff[j];
				ad.im += wj[j].im * this.section[i].bcoeff[j];
			}

			r1 = Math.sqrt(an.re*an.re+an.im*an.im);
			r2 = Math.sqrt(ad.re*ad.re+ad.im*ad.im);

			if ((r1==0.0) || (r2==0.0)) {
				t.re = 0.0;
				t.im = 0.0;
			}
			else {
				t.re = (r1/r2) * ((an.im/r1)*(ad.im/r2) + (an.re/r1)*(ad.re/r2));
				t.im = (r1/r2) * ((an.im/r1)*(ad.re/r2) - (an.re/r1)*(ad.im/r2));
			}

			gain = h.re*t.re - h.im*t.im;
			h.im = h.re*t.im + h.im*t.re;
			h.re = gain;
		}
		gain = Math.sqrt(h.re*h.re + h.im*h.im);
		return(gain);
	}

};
