/*
 * FFT routine based on one originally written by N. Brenner of Lincoln
 * Laboratories.  
 */
 
function RealFFT()
{
	var self=this;
	var me=this.constructor.prototype;

	me.four1 = function(data,nn)
	{
		var	n,mmax,m,j,istep,i;
		var	wtemp,wr,wpr,wpi,wi,theta; /* double precision for the */
                                          		/* trigonometric recurrences. */
		var tempr,tempi;

		n=nn << 1;
		j=1;

		/* Bit-Reversal routine */
		for (i = 1; i < n; i += 2) {
			if (j > i) {
				/* exchange the two complex numbers */
				tempr=data[j]; data[j]=data[i]; data[i]=tempr;
				tempi=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempi;
			}
			m = n >> 1;
			while (m >= 2 && j > m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}                 

		/* Danielson-Lanczos loops */
		mmax = 2;
		while (n > mmax) {      
			/* outer loop executed ln(nn) times */
			istep = 2*mmax; 
                  
			/* initialise for the trigonometric recurrence */
			theta = 2*Math.PI/(mmax); 
			wtemp = Math.sin(0.5*theta);
			wpr   = -2.0*wtemp*wtemp;
			wpi   = Math.sin(theta);
			wr    = 1.0;
			wi    = 0.0;
              
			/* two nested inner loops */
			for (m = 1; m < mmax; m += 2) {
				for (i = m; i <= n; i += istep) {   
					/* Danielson-Lanczos formula: */
					j          = i+mmax;
					tempr      = (wr*data[j]-wi*data[j+1]);
					tempi      = (wr*data[j+1]+wi*data[j]);
					data[j]    = data[i]-tempr;
					data[j+1]  = data[i+1]-tempi;
					data[i]   += tempr;
					data[i+1] += tempi;
				}                                 

				/* trigonometric recurrence */
				wr = (wtemp=wr)*wpr-wi*wpi+wr;
				wi = wi*wpr+wtemp*wpi+wi;
			}
			mmax = istep;
		}
	}

/*
 * This routine takes a single real function and calculates its DFT
 *
 * On Entry: 'sig[0..2n-1]' contains the real function of length '2n',
 * On Exit:  'mag[0..n-1]' contains the magnitudes of the lower half of the DFT.
 */

	me.RealFFTPower = function(sig,len)
	{
		var	i,i1,i2,i3,i4,n2p3;
		var c1=0.5,c2,h1r,h1i,h2r,h2i;
		var	wr,wi,wpr,wpi,wtemp,theta;
		var data, mag;
		var n=2;
		
		/* force length to be power of two */
		while (n < len) n=n*2;
		
		/* copy data into temp buffer */
		data = new Array(n+1);
		for (i=1;i<=len;i++) data[i] = sig[i-1];
		for (;i<=n;i++) data[i]=0;
		n = n >> 1;
		mag = new Array(n);

		/* Initialise the recurrence. */
		theta = Math.PI/n;
		c2 = -0.5; 
		me.four1(data,n);

		wtemp = Math.sin(0.5*theta);
		wpr   = -2.0*wtemp*wtemp;
		wpi   = Math.sin(theta);
		wr    = 1.0+wpr;
		wi    = wpi;
		n2p3  = 2*n+3;  

		for (i = 2; i <= n/2; i++) {        /* Case i==1 done below */

			/* Separate transforms from the data */
			i4  = 1+(i3=n2p3-(i2=1+(i1=i+i-1)));
			h1r = c1*(data[i1]+data[i3]);
			h1i = c1*(data[i2]-data[i4]);
			h2r = -c2*(data[i2]+data[i4]);
			h2i = c2*(data[i1]-data[i3]);

			/* Recombine to form true transform of original data */
			data[i1] = (h1r+wr*h2r-wi*h2i);
			data[i2] = (h1i+wr*h2i+wi*h2r);
			data[i3] = (h1r-wr*h2r+wi*h2i);
			data[i4] = (-h1i+wr*h2i+wi*h2r);

			/* The recurrence Relations */
			wr = (wtemp = wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}    

		/* Squeeze the first & last data into the original array */ 
		data[2]=0;
		data[1] = ((h1r=data[1])+data[2]);
		data[2] = (h1r-data[2]);

		for (i=0;i<n;i++) {
			mag[i] = Math.sqrt(data[2*i+1]*data[2*i+1]+data[2*i+2]*data[2*i+2])/n;
		}
		
		return mag;
	} 
}
