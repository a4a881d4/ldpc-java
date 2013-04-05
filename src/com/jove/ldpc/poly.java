package com.jove.ldpc;

import java.io.*;

/**
 * 
 */

/**
 * @author zhaom
 *
 */


public class poly {
	int polys[];
	int length;
	public static PrintWriter out;
	public static int M;
	public static int GP;
	public static int P;
	static int[] index;
	static int[] power;
	
	int gmul( int i, int j )
	{
		if( i==0 || j==0 )
			return 0;
		int res=index[i]+index[j];
		res%=P;
		return power[res];
	}

	int gdiv( int i, int j )
	{
		if( i==0 )
			return 0;
		int res=index[i]-index[j];
		if( res<0 )
			res+=P;
		return power[res];
	}

	int gadd( int a, int b)
	{
		return(a^b)&P;
	}

	public static void init(int gp)
	{
		poly.GP = gp;
		int i;
		i = gp;
		poly.M = 0;
		for(;i!=1&&i!=0;i>>=1) poly.M++;
		poly.P = (1<<poly.M)-1;
		poly.index=new int[P+1];
		poly.power=new int[P+1];
		
		int a=1;
		{
			for( i=0;i<P;i++ )
			{
				power[i]=a;
				index[a]=i;
				a<<=1;
				if( ((a>>M)&1)==1 )
					a^=GP;
			}
		}
		index[0]=P;
		power[P]=0;
	}

	poly(int MaxLen)
	{
		length = 0;
		polys = new int[MaxLen];
	}
	public poly()
	{
		length = 0;
		polys = new int[1024];
	}
	
	public poly copyPoly( poly a)
	{
		poly b = new poly();
		b.length = a.length;
		int i;
		for( i=0;i<a.length;i++ )
			b.polys[i]=a.polys[i];
		return b;
	}
	public void copyPoly( poly a, poly b)
	{
		b.length = a.length;
		int i;
		for( i=0;i<a.length;i++ )
			b.polys[i]=a.polys[i];
	}
	public void print( String name )
	{
		int i;
		int j=0;
		out.println("Polynomail "+name+":");
		for( i=length-1;i>0;i-- )
		{
			if(polys[i]!=0)
			{
				out.println(polys[i]+"X**"+i+"+");
				j++;
			}
		}
		if( polys[0]!=0 && length!=0 )
			out.println(polys[0]);
	}
	
	public poly xk( int k )
	{
		poly a = new poly();
		a.length=k+1;
		for( int i=0;i<k;i++ )
		{
			a.polys[i]=0;
		}
		a.polys[k]=1;
		return a;
	}
	
	public poly pdiv( poly a, poly b )
	{

		poly d = new poly();
		int i,j;
		d.length=a.length+1-b.length;
		//assert( b.polys[b.length-1]==1 );
		for( i=a.length-1;i>=b.length-1;i-- )
		{
			d.polys[i-b.length+1]=a.polys[i];
			for( j=0;j<b.length-1;j++ )
				a.polys[i-j-1]^=gmul(a.polys[i],b.polys[b.length-2-j]);
		}
		a.length=b.length-1;
		for( i=a.length-1;i>=0;i-- )
		{
			if( a.polys[i]!=0 )
			{
				a.length=i+1;
				break;	
			}
		}
		if( i==-1 )
			a.length=0;
		return d;
	}
	public poly padd( poly a, poly b )
	{
		poly d;
		int i;
		if( a.length>b.length )	
		{
			d = copyPoly(a);
			for( i=0;i<b.length;i++ )
				d.polys[i]^=b.polys[i];	
		}
		else
		{
			d = copyPoly(b);
			for( i=0;i<a.length;i++ )
				d.polys[i]^=a.polys[i];	
		}
		for( i=d.length-1;i>=0;i-- )
		{
			if( d.polys[i]!=0 )
			{
				d.length=i+1;
				break;	
			}
		}
		if( i==-1 )
			d.length=0;
		return d;
	}
	
	public poly pmul( poly a, poly b )
	{

		int i,j;
		poly d = new poly();
		d.length=a.length+b.length-1;
		for( i=0;i<d.length;i++)
			d.polys[i]=0;
		for( i=0;i<b.length;i++ )
		{
			for( j=0;j<a.length;j++ )
				d.polys[i+j]^=(gmul(a.polys[j],b.polys[i]));
		}
		return d;
	}

	
	public int rdiv( poly a, poly b, poly m, poly n, int debug )
	{
		if( debug==1 )
		{
			a.print( "a");
			b.print( "a");
				
		}
		int i,j;
		poly[] dl = new poly[127];
		poly g1 = new poly();
		poly g2;
		poly g3;
		poly pg1,pg2,pg3;
		poly mk0,mk1,mk2;
		
		g3 = copyPoly(b);
		g2 = copyPoly(a);
		pg1=g1;
		pg2=g2;
		pg3=g3;
		i=0;
		poly t;
		
		while( pg3.length != 1 )	
		{
			t=pg1;
			pg1=pg2;
			pg2=pg3;
			pg3=t;
			if( debug==1 )
			{
				pg1.print( "div1");
				pg2.print("div2");
				
			}
			dl[i]=pdiv(pg1,pg2);
			
			if( debug==1 )
			{
				dl[i].print("d["+i+"]");
				pg1.print("rem");
			}
			pg3 = copyPoly(pg1);
			i++;
			if( pg3.length==0 )
			{
				copyPoly(pg2,m);
				return -1;
			}
		}
		
		mk2=g3;
		mk1=g2;
		mk0=g1;
		
		mk2.length=1;
		mk2.polys[0]=1;
		
		mk1 = copyPoly( dl[i-1] );
		copyPoly( mk2,n );
		copyPoly( mk1,m );
		
		for( j=i-2;j>=0;j-- )
		{
			if( debug==1 )
			{
		
				mk2.print("mk2");
				mk1.print("mk1");
				dl[j].print("dl["+j+"]");
			}
			mk0 = pmul( mk1,dl[j]);
			if( debug==1 )
			{
				mk0.print( "before add mk0");
			}
			mk0 = padd( mk0, mk2 );
			if( debug==1 )
			{
				mk0.print( "mk0");
			}
			copyPoly( mk1,n );
			copyPoly( mk0,m );
			t=mk2;
			mk2=mk1;
			mk1=mk0;
			mk0=t;
		}
		
		if( debug==1 )
		{	
			m.print( "m");
			n.print( "m");
		}
		return 0;	
			
	}
	
	static poly pinv(poly a, int k)
	{
		poly g = a.xk(k);;
		g.polys[0]=1;
		poly m = new poly();
		poly n = new poly();
		a.rdiv(a,g,m,n,0);
		poly b = a.padd(a.pmul(a, n),a.pmul(g, m));
//		b.print("b:");
		if(b.length!=1)
			return null;
		return n;
	}
	
	public int findp( poly[] dpp, int p )
	{
		poly a,b,c,d;
		int i,l,k;
		
		a=xk(p);
		a.polys[0]=1;

		i=3;
		l=0;
	
		while(i<0x8000000)
		{
			b = new poly();
			k=i;
			b.length=0;
			for(;k!=0;k>>=1)
			{
				b.polys[b.length]=(char)(k&1);
				b.length++;
			}
			for(k=0;k<l;k++)
			{
				d = copyPoly(b);
				c = pdiv(d,dpp[k]);
				if( c.length==0 )
					break;
			}
			while( k==l )
			{
				d = copyPoly(a);
				c = pdiv(d,b);
				if( d.length!=0 )
					break;
				dpp[l] = copyPoly(b);
				l++;
				copyPoly(c,a);
			}
			if( a.length<b.length )
				break;
			i++;
		}
		return l;	
	}
	
	public int findp( Integer[] ipp, int p )
	{
		int i,j,k;
		poly[] dp = new poly[127];
		int num = findp( dp, p );
		for( i=0;i<num;i++ )
		{
			k=0;
			for(j=0;j<dp[i].length;j++)
				if(dp[i].polys[j]==1)
					k|=(1<<j);
			ipp[i]=k;
		}
		return num;
	}
	
	void genGPoly( int t )
	{
		poly b = new poly();
		poly c = new poly();
		c.length=2;
		c.polys[0]=power[0];
		c.polys[1]=1;
		int i;
		for( i=1;i<=t;i++ )
		{
			b.length=2;
			b.polys[0]=power[i];	
			b.polys[1]=1;
			c = pmul( c, b );
		}
		length = c.length;
		for( i=0;i<c.length;i++ )
		{
			polys[i]=c.polys[i];
		}
	}

	void int2Poly( int[] a, int length )
	{
		int i;
		this.length = length; 
		for(i=0;i<length;i++)
		{
			polys[i]=(a[i/32]>>(i%32))&1;
		}
	}
	int poly2Int(int[] buf)
	{
		int i;
		assert(buf.length>=(length+31)/32);
		for(i=0;i<buf.length;i++) buf[i]=0;
		for(i=0;i<length;i++)
		{
			if( polys[i]!=0 )
				buf[i/32]|=1<<(i%32);
		}
		return length;
	}
	public int[][] QCRSTab()
	{
		int [][] ret = new int[poly.P][];
		genGPoly(poly.P-3);
		int i,j;
		for( i=0;i<poly.P;i++ )
		{
			poly d = new poly();
			d.length = 2;
			d.polys[1]=poly.power[0];	
			d.polys[0]=poly.power[i];
			poly b = pmul( d, this );
			ret[i] = new int[b.length];
			for( j=0;j<b.length;j++ )
			{
				ret[i][j] = poly.index[b.polys[j]];
			}
		}
		return ret;
	}
	public static void main(String[] args)
	{
		int i,j;
		int k;
		int p;
		poly.out = new PrintWriter(System.out,true);
		poly.init(0x3);
		p = 127;
		
		poly[] dp = new poly[127];
		poly[] mp = new poly[127];
		poly m = new poly();
		poly n = new poly();
		poly a,b,c,e;
		
		int num = m.findp( dp, p );	
		
		for( i=0;i<num;i++ )
		{
			a = m.xk(0);
			for( j=0;j<num;j++ )
			{
				if( i!=j )
				{	
					a = a.pmul(a,dp[j]);
				}
			}
			m.length=0;
			n.length=0;
			a.rdiv(a,dp[i],m,n,1);
			e = a.pmul(n,a);
			b = a.pmul(m,dp[i]);
			c = a.padd(e,b);
			if( c.length!=1 && c.polys[0]!=1 )
			{
				c.print("Error c");
			}
			mp[i] = a.copyPoly(n);
		}
		for( i=0;i<num;i++ )
		{
			k=0;
			for(j=0;j<dp[i].length;j++)
				if(dp[i].polys[j]==1)
					k|=(1<<j);
			Integer Ik = new Integer(k);
			poly.out.printf("0x%x, ",Ik);
			k=0;
			for(j=0;j<mp[i].length;j++)
				if(mp[i].polys[j]==1)
					k|=(1<<j);
			Ik = new Integer(k);
			poly.out.printf("0x%x, \n",Ik);
		}
		Integer[] ipp = new Integer[128];
		num = m.findp( ipp, 127);
		for( i=0;i<num;i++ )
			poly.out.printf("0x%x ", ipp[i]);
		poly.out.printf("\n");
		
		poly.out.println("Try inv");
		
		poly g = m.xk(136);
		g.polys[0]=1;
		for( i=0;i<num;i++ )
		{
			poly invdp = poly.pinv(dp[i], 136);
			if( invdp == null )
				poly.out.println("dp["+i+"] can not inv");
			else
			{
				poly product = invdp.pmul(dp[i], invdp);
				product.pdiv(product,g);
				product.print("Product :");
			}
		}

		poly.init(0x43);
		a = new poly();
		int[][] rs = a.QCRSTab();
		
		poly.out.println(rs.length);
		
		
	}
}
