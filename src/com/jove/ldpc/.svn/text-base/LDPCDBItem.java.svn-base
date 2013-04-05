package com.jove.ldpc;

import java.io.PrintWriter;

public class LDPCDBItem {
	Degree deg;
	public String comments;
	String Url;
	public double rate;
	public double[] EbN0Range;
	double[] Pb;
	public int[][] code;
	LDPCDBItem( Degree deg, String comm, String url, double[] Eb, double rate )
	{
		this.deg = deg;
		comments = comm;
		Url = url;
		this.rate = rate;
		EbN0Range = Eb;
	}
	void setEb( double[] Eb )
	{
		this.EbN0Range = Eb;
	}
	public int[][] getCode( int N, int bl, int poly, PrintWriter print )
	{
		int i,j;
		int M = (int)Math.round((double)N * (1.-rate));
		QCRS rs = new QCRS( deg, M, N, bl, poly );
		int pos=0;
		rs.genH(0,0);
		while( rs.calcPhyInv(0)==false )
		{
			rs.swap(rs.H,pos,(rs.N-rs.M));
			System.out.println("Swap "+pos+" "+(rs.N-rs.M));
			pos++;
			if( pos> rs.N-M-1 )
				break;
		};
		if( pos<=rs.N-M-1 )
		{
			if( print != null )
			{
				print.printf("%d %d %d<br>",rs.N,rs.M,rs.bl);
				print.printf("<br>int[][] rs = {<br>");
				for( i=0;i<rs.M;i++ )
				{
					print.printf("{ ");
					for( j=0;j<rs.N;j++ )
						print.printf("%4d, ",rs.H[i][j]);
					print.printf("},<br>");
				}
				print.printf("};<br>");
				
			}
			code = rs.H;
			return rs.H; 
		}
		else
			return null;
	}
	public float[] sigma()
	{
		float[] ret = new float[EbN0Range.length];
		for(int i=0;i<EbN0Range.length;i++ )
		{
			ret[i] = 0.5f/(float)(Math.pow(10.,EbN0Range[i]/10.)*rate);
			ret[i] = (float)Math.sqrt(ret[i]);
		}
		return ret;
	}
}
 