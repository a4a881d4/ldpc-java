package com.jove.ldpc;

public class AWGN {
	
	float var;
	
	public AWGN( float var )
	{
		this.var = (float)var;
		c=0;
	}
	static float calcN(float[] a)
	{
		float n = 0f;
		for( int i=0;i<a.length;i++ )
		{
			if( a[i]>0f )
				n+= (a[i]-1f)*(a[i]-1f);
			else
				n+= (a[i]+1f)*(a[i]+1f);
		}
		n/=(float)a.length;
		n = (float)Math.sqrt((double)n);
		return n;
	}
	int c;
	float mem;
	private float randn( )
	{
		float now;
		if( c==0 )
		{
			c=1;
			float a = (float)Math.random();
			float b = (float)Math.random();
			now = (float)(Math.cos(2.0*Math.PI*a) * Math.sqrt(-2.0*Math.log(b)));
			mem = (float)(Math.sin(2.0*Math.PI*a) * Math.sqrt(-2.0*Math.log(b)));
			return now*var;
		}
		else
		{
			c=0;
			return mem*var;
		}
	}
	private float channel(int i)
	{
		float n = randn();
		float ret = (i==0)? (float)1.+n:(float)-1.+n;
		return ret;
	}
	public float[] channel( byte[] in)
	{
		float[] ret = new float[in.length*8];
		for( int i=0;i<in.length*8;i++ )
		{
			ret[i] = channel((in[i/8]>>(i%8))&1);
		}
		return ret;
	}
	int channel( float[] out, int k, char[] in)
	{
		for( int i=0;i<in.length*8;i++ )
		{
			out[i+k] = channel((in[i/8]>>(i%8))&1);
		}
		return in.length*8;
	}
}
