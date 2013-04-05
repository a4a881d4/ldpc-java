package com.jove.ldpc;

import com.jove.ldpc.LDPCDBItem;
import java.io.PrintWriter;

public class Dde {
	/*
  	FILENAME: dderegular.c
  	AUTHOR: Tadashi Wadayama

	HOW TO MAKE:
	gcc -O2 -o dderegular dderegular.c -lm -lfftw3

	The program performs density evolution for (quantized) BP over AWGN channel.
	The program is based on the idea of "discretized DE" by S.Y.Chung.
	A (j,k)-regular LDPC ensemble is assumed.
	This program uses  FFTW FFT library (see http://www.fftw.org/).

	Note:

	array      <-> pmf
	---------------------------------------------
	a[0]       <-> a_{-(ASIZE-1)/2}  * DELTA
	a[1]       <-> a_{-(ASIZE-1)/2+1}* DELTA
	a[2]
	... 
	a[(ASIZE-1)/2] <-> a_0 *  DELTA
	... 
	a[ASIZE-1] <-> a_{+(ASIZE-1)/2} * DELTA
	---------------------------------------------

	zero padding
	---------------------------------------------
	a_LB = a_{LB+1} = ... = a_{LB2-1} = 0
	a_UB = a_{UB-1} = ... = a_{UB2+1} = 0
	---------------------------------------------
	LB  = -(ASIZE-1)/2
	UB  =  (ASIZE-1)/2

	LB2 = -(ASIZE-1)/4
	UB2 =  (ASIZE-1)/4
	---------------------------------------------

	pdf(discrete) <-> pdf (real value)
	---------------------------------------------
	a_i       <-> [DELTA i- DELTA/2, DELTA i + DELTA/2]
	(DELTA: step size)
	---------------------------------------------

  	Version 1.1

	Copyright (C) Tadashi Wadayama
	June 28, 2003: start
	Aug.6: fixed: definition of pi, quantize_pdf(initialization bug), apply_R(normalization bug)
	Aug.7: fixed: error_prob

	TODO:
	improve efficiency of Rpower()
	 */

	static class pdf {
		double x;
		static pdf[] alloc(int size)
		{
			pdf[] r = new pdf[size];
			for( int i=0;i<size;i++)
				r[i] = new pdf();
			return r;
		}
		static void clean( pdf[] a )
		{
			for( int i=0;i<a.length;i++ )
				a[i].x = 0.0;
		}
		static void print_pdf(pdf[] a,int len, String name )
		{
			int i;
			int mid = (Dde.ASIZE-1)/2;
			len = (len-1)/2;
			System.out.println("PDF["+name+"]:");
			for (i = mid-len; i <= mid+len; i++) 
				System.out.printf("%16.12f %16.12f\n",DELTA*(i-mid),a[i].x);
		}

	}

	static int LB,UB,LB2,UB2;
	static double UV, LV;

	static class fft_plan {
		Complex[] in;
		Complex[] out;
		int inv;
		int N;
		static int count;
		fft_plan(int N,Complex[] in,Complex[] out, int inv )
		{
			this.in = in;
			this.out = out;
			this.N =N;
			this.inv = inv;
			count = 0;
		}
		void doit()
		{
			if(inv==0)
				out = FFT.fft(in);
			else
				out = FFT.ifft(in);
			count++;
		}
	}

	fft_plan p1,p2,p3;

	Complex[] _a, _b, _c, _A, _B, _C;
	pdf[]  tmp1, tmp2, tmp3, tmp4;
	int[][] R_tbl;

	static int ASIZE;			/* array size */
	static double DELTA;			/* step size of quantization */


	double channel_var;

	void clip_pdf(pdf[] x)
	/* PDF clipping function */
	{
		int i;
		double sum;
		sum = 0.0;
		for (i = LB+UB; i <= UB+LB2-1; i++) {
			sum += x[i].x;
			x[i].x = 0.0;
		}
		x[LB2+UB].x +=sum;
		sum = 0.0;
		for (i = UB+UB2+1; i <= ASIZE-1; i++) {
			sum += x[i].x;
			x[i].x = 0.0;
		}
		x[UB2+UB].x +=sum;
	}

	double llr_pdf(double x)
	/* PDF of LLR (AWGN channel) */
	{
		double var = 4.0/channel_var;
		double mean = 2.0/channel_var;
		return (1.0/Math.sqrt(2.0*Math.PI*var)) * Math.exp(-(x-mean)*(x-mean)/(2.0*var));
	}


	void quantize_pdf(pdf[] a)
	/* The function makes quatized version of a given PDF function.  */
	{
		int i;
		double l,r;
		double vl,vr;


		for (i = LB2; i <= UB2; i++) {
			l = DELTA * i - DELTA/2.0;
			r = DELTA * i + DELTA/2.0;
			vl = llr_pdf(l);
			vr = llr_pdf(r);
			a[i+UB].x = (vl+vr)/2.0;
		}
		for (i = LB; i <= LB2-1; i++) a[i+UB].x = 0.0;
		for (i = UB2+1; i <= UB; i++) a[i+UB].x = 0.0;
	}



	double integrate_pdf(pdf[] a)
	/* The function integrates a given PDF. */
	{
		int i;
		double sum;
		sum = 0.0;
		for (i = LB; i <= UB; i++) {
			sum += a[i+UB].x;
		}
		return sum * DELTA;
	}

	double error_prob(pdf[] a)
	{
		int i;
		double sum;
		sum = 0.0;
		for (i = 0;i <= UB-1; i++) {
			sum +=  a[i].x;
		}
		sum += a[UB].x/2.0;
		return sum * DELTA;
	}

	void normalize_pdf(pdf[] r)
	{
		int i;
		double sum;

		sum = 0.0;			/* normalization */

		for (i = 0; i <= ASIZE-1; i++) sum += r[i].x*DELTA;
		for (i = 0; i <= ASIZE-1; i++) r[i].x/= sum;
	}

	void conv(pdf[] a, pdf[] b, pdf[] r)
	/* convolution of two PDFs */
	{
		int i;
		//  double sum;

		for (i = 0; i < ASIZE; i++) 
		{
			_a[i] = new Complex(a[i].x,0.);
			_b[i] = new Complex(b[i].x,0.0);
		}
		_a[ASIZE] = new Complex(0.0,0.0);
		_b[ASIZE] = new Complex(0.0,0.0);

		p1.doit();
		p2.doit();
		for (i = 0; i < ASIZE+1; i++) 
		{
			p3.in[i] = p1.out[i].times(p2.out[i]);
		}
		p3.doit();
		for (i = 0; i <= ASIZE-1; i++) 
		{
			r[i].x = p3.out[((ASIZE-1)/2+i)%(ASIZE+1)].re() * DELTA;
		}
		clip_pdf(r);
	}

	void delta_func(pdf[] a)
	/* delta function */
	{
		int i;
		for (i = 0; i <= ASIZE-1; i++) {
			a[i].x = 0.0;
		}
		a[UB].x = 1.0/DELTA;
	}

	void naive_convpower(pdf[] a, pdf[] b, int n)
	{
		int i;
		for (i = 0; i <= ASIZE-1; i++) {
			b[i].x = a[i].x;
		}
		for (i = 1; i <= n-1; i++) {
			conv(a,b,b);
		}
	}

	void convpower(pdf[] a, pdf[] b, int n)
	{
		int i;
		int w;
		int bit;

		if (n < 6) {
			naive_convpower(a,b,n);
			return;
		}
		for (i = 0; i <= ASIZE-1; i++) tmp1[i].x = a[i].x;
		delta_func(b);
		w = 0;
		i = 0;
		while( true ) {
			bit = (n >> i) & 1;
			if (bit == 1) conv(b,tmp1,b);
			if (bit == 1) w += (1<<i);
			if (w == n) break;
			i++;
			conv(tmp1,tmp1,tmp1);
		}
	}
	double atanh(double x)
	{
		return 0.5*Math.log((1+x)/(1-x));
	}

	double R(double a, double b)
	{
		return 2.0 * atanh(Math.tanh(a/2.0)*Math.tanh(b/2.0));
	}

	void init()
	/* initialization */
	{
		int i,j;
		int idx;
		double a,b,v;

		ASIZE = (1 << 12)-1;
		DELTA = 50.0/ASIZE;

		LB = -(ASIZE-1)/2;		
		UB = (ASIZE-1)/2;
		LB2 = -(ASIZE-1)/4;		
		UB2 = (ASIZE-1)/4;
		UV = UB * DELTA + DELTA/2.0;
		LV = LB * DELTA - DELTA/2.0;

//		System.out.printf("#ASIZE= %d\n",ASIZE);
//		System.out.printf("#LB  = %d\n",LB);
//		System.out.printf("#UB  = %d\n",UB);
//		System.out.printf("#LB2 = %d\n",LB2);
//		System.out.printf("#UB2 = %d\n",UB2);
//		System.out.printf("#UV  = %f\n",UV);
//		System.out.printf("#LV  = %f\n",LV);

		_a = new Complex[ASIZE+1];
		_b = new Complex[ASIZE+1];
		_c = new Complex[ASIZE+1];
		_A = new Complex[ASIZE+1];
		_B = new Complex[ASIZE+1];
		_C = new Complex[ASIZE+1];

		tmp1 = pdf.alloc(ASIZE);
		tmp2 = pdf.alloc(ASIZE);
		tmp3 = pdf.alloc(ASIZE);
		tmp4 = pdf.alloc(ASIZE);

		R_tbl = new int[ASIZE][];
		for (i = 0; i < ASIZE; i++) 
		{
			R_tbl[i] = new int[ASIZE];
		}

		for (i = LB2; i <= UB2; i++) {
			for (j = LB2; j <= UB2; j++) {
				a = DELTA * i;
				b = DELTA * j;
				v = R(a,b);
				idx = 0;
				if (v >= DELTA/2.0) 
					idx = (int)Math.floor(v/DELTA + 0.5);
				if (v <= -DELTA/2.0) 
					idx = (int)Math.ceil(v/DELTA - 0.5);
				if ((v < DELTA/2.0) && (v > -DELTA/2.0)) 
					idx = 0;
				if (v > UB2*DELTA) 
					idx = UB2;
				if (v < LB2*DELTA) 
					idx = LB2;
				R_tbl[i+UB][j+UB] = idx;
			}
		}
		p1 = new fft_plan(ASIZE+1, _a, _A, 0 );  
		p2 = new fft_plan(ASIZE+1, _b, _B, 0 );  
		p3 = new fft_plan(ASIZE+1, _C, _c, 1 );  
	}

	void apply_R(pdf[] a, pdf[] b, pdf[] c)
	{
		int i,j;
		for (i = 0; i <= ASIZE-1; i++) 
			tmp1[i].x = 0.0;
		for (i = LB2+UB; i <= UB2+UB; i++) 
		{
			for (j = LB2+UB; j <= UB2+UB; j++) 
			{
				tmp1[R_tbl[i][j]+UB ].x += a[i].x * b[j].x;
			}
		}
		for (i = 0; i <= ASIZE-1; i++) 
			c[i].x = tmp1[i].x * DELTA;
		
	}

	void Rpower(pdf[] a, pdf[] b, int n)
	{
		int i;
		for (i = 0; i <= ASIZE-1; i++) 
			b[i].x = a[i].x;

		for (i = 1; i <= n-1; i++)
		{	
			//pdf.print_pdf(b, 11, "Rpower ["+i+"]");
			apply_R(a,b,b);
		}
			
	}

	void Rpower_ef(pdf[] a,pdf[] b, int n)
	{
		int i,j;
		int w;
		int bit;

		if (n < 6) {
			Rpower(a,b,n);
			return;
		}
		for (i = 0; i <= ASIZE-1; i++) tmp1[i] = a[i];
		// delta_func(b);
		w = 0;
		i = 0;
		while( true ) {
			bit = (n >> i) & 1;
			if (bit == 1)
			{
				if (i==0)
				{
					for (j=0;j<ASIZE-1;j++)
						b[j]=a[j];
				}
				else
				{
					apply_R(b,tmp1,b);
				}
			}
			if (bit == 1) w += (1<<i);
			if (w == n) break;
			i++;
			apply_R(tmp1,tmp1,tmp1);
		}	
	}

	void Add_fun(pdf[] f0,pdf[] f1,double fraction)
	{
		int i;
		for(i=0;i<ASIZE;i++)
		{
			f0[i].x+=fraction*f1[i].x;
		}
	}
class EvolutionRet {
	int iteration;
	int stop;
	double pb;
	EvolutionRet(int iteration, int stop, double pb )
	{
		this.iteration = iteration;
		this.stop = stop;
		this.pb = pb;
	}
}
EvolutionRet evolution(Degree lamda, Degree beta, double var, int iteration, int print  )
	{
		channel_var = var*var;
		int i,j;
		double pb;
		double pbMem;
		pdf[] P_lambda, P_c2m, P_m2c, f0, f1, tmpA, tmpB, tmpC, tmpD;
		P_lambda = pdf.alloc(ASIZE);
		P_c2m = pdf.alloc(ASIZE);
		P_m2c = pdf.alloc(ASIZE);
		f0 = pdf.alloc(ASIZE);
		f1= pdf.alloc(ASIZE);
		tmpA = pdf.alloc(ASIZE);
		tmpB = pdf.alloc(ASIZE);
		tmpC = pdf.alloc(ASIZE);
		tmpD = pdf.alloc(ASIZE);

		/* initialization of LLR density */
		/* initialization of message(c->m) density */
		quantize_pdf(P_lambda);
		normalize_pdf(P_lambda);
		delta_func(P_c2m);
		j=0;
		pb=1.;
		pbMem = pb;
		while(j<iteration) 
		{

			/* variable node operation */
			tmpC = pdf.alloc(ASIZE);
			for (i=0;i<lamda.tab.length;i++)
			{
				convpower(P_c2m, tmpA, lamda.tab[i].deg-1);
				normalize_pdf(tmpA);
				Add_fun(tmpC,tmpA,lamda.tab[i].degFrac);
			}
			conv(P_lambda, tmpC, P_m2c);
			normalize_pdf(P_m2c);

			conv(P_m2c,P_c2m,tmpC);
			pb = error_prob(tmpC);
			if( ((print>>0)&1) == 1 )
				System.out.printf("# %d %16.12f\n",j,pb);
			if (Math.abs(pb) < 1e-6) 
				return new EvolutionRet(j,0,pb); 
			if(Math.abs(pb-pbMem)<1e-10)
				return new EvolutionRet(j,-1,pb); 
			pbMem = pb;	
			/* check node operation */

			pdf.clean(P_c2m);
			for (i=0;i<beta.tab.length;i++)
			{
				Rpower(P_m2c,tmpC,beta.tab[i].deg-1);
				normalize_pdf(tmpC);
				Add_fun(P_c2m,tmpC,beta.tab[i].degFrac);
			}
			j++;
		}
		return new EvolutionRet(j,0,pb); 
	}
	Dde() 
	{
		init();
	}
	EvolutionRet ret;
	public Dde( int idx, double EbN0,PrintWriter out )
	{
		init();
		LDPCDB adb = new LDPCDB();
		LDPCDBItem x = adb.db.get(idx);
		double rate = x.rate;
		Degree beta = x.deg.v2c(rate).node2edge();
		Degree lamda = x.deg.node2edge();
		out.println("Code rate: "+Degree.rateE(lamda, beta));
		double sigma = Math.sqrt((1./rate)*0.5/Math.pow(10., EbN0/10.));
		ret = evolution(lamda, beta, sigma, 1000, 0 );
		out.println(x.comments);
		out.println("Eb/N0 = "+EbN0+" Noise "+sigma+" After Iteration "+ret.iteration+" the Pb: "+ret.pb);
	}
	public static void main(String[] args)
	{
		
		LDPCDB adb = new LDPCDB();
		for( LDPCDBItem x : adb.db )
		{
			double rate = x.rate;
			Degree beta = x.deg.v2c(rate).node2edge();
			Degree lamda = x.deg.node2edge();
			System.out.println("Code rate: "+Degree.rateE(lamda, beta));
			Dde ttu = new Dde();
			for( double EbN0 : x.EbN0Range )
			{
				double sigma = Math.sqrt((1./rate)*0.5/Math.pow(10., EbN0/10.));
				EvolutionRet ret = ttu.evolution(lamda, beta, sigma, 1000, 0 );
				System.out.println(x.comments);
				System.out.println("Eb/N0 = "+EbN0+" Noise "+sigma+" After Iteration "+ret.iteration+" the Pb: "+ret.pb);
				if( ret.stop==-1 )
					break;
			}
		}
	}
}

