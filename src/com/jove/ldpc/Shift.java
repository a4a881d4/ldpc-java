package com.jove.ldpc;

public class Shift {

	/**
	 * @param args
	 */
	int[] buf;
	static int length;
	static int size;
	static int mask;
	
	static void init()
	{
		Shift.size = (length+31)/32;
		int left = size*32 -length;
		Shift.mask = shiftMy((0xffffffff),left); 
	}
	
	static int shiftMy( int a, int s)
	{
		if( s>=1 )
		{
			a = a>>1;
			a &= 0x7fffffff;
			a = a>>(s-1);
		}
		return a;
	}
	
	Shift()
	{
		buf = new int[size];
	}
	void clean()
	{
		for(int i=0;i<size;i++)
			buf[i]=0;
	}
	static void clean( Shift[] a)
	{
		for( int i=0;i<a.length;i++ )
			a[i].clean();
	}
	void dump(Shift a)
	{
		for(int i=0;i<size;i++)
			buf[i]=a.buf[i];
	}
	static void cpy(Shift[] from, Shift[] to)
	{
		for( int i=0;i<to.length;i++ )
			to[i].dump(from[i]);
	}
	int weight()
	{
		int i,j;
		i=0;
		for( j=0;j<length;j++ )
			if( ((buf[j/32]>>(j%32))&1) == 1 )
				i++;
		return i;
	}
	void mac(Shift a,int s)
	{
		Shift b = a.shiftH2L(s);
		for(int i=0;i<size;i++)
			buf[i]^=b.buf[i];
	}
	void setXk( int k )
	{
		this.zero();
		buf[k/32]=1<<(k%32);
	}
	void zero()
	{
		for( int i=0;i<size;i++)
			buf[i]=0;
	}
	static void zero(Shift[] a)
	{
		for( int i=0;i<a.length;i++ )
			a[i].zero();
	}
	static Shift[] alloc(int k)
	{
		Shift[] ret = new Shift[k];
		for( int i=0;i<k;i++ )
			ret[i]= new Shift();
		Shift.zero(ret);
		return ret;
	}
	static Shift[][] alloc(int[][] a)
	{
		int i,j;
		Shift[][] ret = new Shift[a.length][];
		for( i=0;i<a.length;i++ )
		{
			ret[i]=Shift.alloc(a[i].length);
			for( j=0;j<a[i].length;j++ )
			{
				if(a[i][j]!=Shift.length)
				{
					ret[i][j].setXk((Shift.length-a[i][j])%Shift.length);
				}
			}
		}
		return ret;
	}
	static Shift reverse(Shift a)
	{
		Shift b = new Shift();
		b.zero();
		for(int i=0;i<length;i++)
		{
			int j = (length-i)%length;
			int c = (a.buf[j/32]>>(j%32))&1;
			b.buf[i/32]|=c<<(i%32);
		}
		return b;
	}
	int isXk()
	{
		int i,j;
		for( i=0;i<length;i++ )
		{
			if( ((buf[i/32]>>(i%32))&1)==1 )
				break;
		}
		if( i==length )
			return -2;
		for( j=i+1;j<length;j++ )
		{
			if( ((buf[j/32]>>(j%32))&1)==1 )
			{
				System.out.printf("%d is one but %d is not zero, \"buf[%d]=%08x, buf[%d]=%08x\"\n",
						i,j,i/32,buf[i/32],j/32,buf[j/32]);
				return -1;
			}
		}
		return i;
	}
	
	int mid( int []a, int begin, int end, int to )
	{
		int r = 0;
		if(begin>end)
		{
			r = mid(a,0,end,to);
			r |= mid(a,begin,length-1,length-1-begin);
			return r;
		}
		int posE = end/32;
		int leftE = end%32;
		int posB = begin/32;
		int leftB = begin%32;
		if(posE==posB)
		{
			r = a[posE]<<(31-leftE);
			r = shiftMy( r , (31-leftE+leftB));
			r = r<<(to-(end-begin));
		}
		else
		{
			r = shiftMy(a[posB],leftB);
			r |= a[posE]<<(to-leftE);
		}
		return r;
	}
	void shiftOneInt(Shift a, int s)
	{
		int size = (length+31)/32;
		int i;
		for( i=0;i<size-1;i++ )
		{
			buf[i]=mid(a.buf,(i*32+s)%length,(i*32+31+s)%length,31);
		}
		buf[size-1]=mid(a.buf,((size-1)*32+s)%length,(length-1+s)%length,(length-1)%32);
		buf[size-1]&=mask;
	}
	Shift shiftOne(int s)
	{
		Shift b = new Shift();
		int i;
		for( i=0;i<size-1;i++ )
		{
			b.buf[i]=mid(buf,(i*32+s)%length,(i*32+31+s)%length,31);
		}
		b.buf[size-1]=mid(buf,((size-1)*32+s)%length,(length-1+s)%length,(length-1)%32);
		b.buf[size-1]&=mask;
		return b;
	}
	
	Shift shiftH2L(int s)
	{
		return shiftOne(s);
	}
	
	Shift shiftL2H(int s)
	{
		return shiftOne(length-s);
	}
	
	Shift shiftOneByExt( int s)
	{
		Shift r = new Shift();
		int i;
		int[] temp = new int[length];
		for( i=0;i<length;i++)
			temp[i] = (buf[i/32]>>(i%32))&1;
		for( i=0;i<size;i++ )
			r.buf[i]=0;
		for( i=0;i<length;i++)
		{
			if( temp[(i+s)%length]== 1 )
				r.buf[i/32]|=1<<(i%32);
		}	
		return r;
	}
	
	boolean Equ( Shift a )
	{
		for( int i=0;i<size;i++ )
			if(a.buf[i]!=buf[i])
				return false;
		return true;
	}
	void print()
	{
		int i;
		int c = 0;
		for( i=0;i<length;i++ )
		{
			c<<=1;
			c|=(buf[i/32]>>(i%32))&1;
			if((i%4)==3)
				System.out.printf("%01x",c&0xf);
		}
	}
	void printH2L()
	{
		int i;
		int c = 0;
		for( i=size-1;i>=0;i-- )
		{
			System.out.printf("%08x",buf[i]);
		}
	}
	static void printArray(Shift[] a, String name, int comma )
	{
		int i;
		
		if((comma&2) == 2 && name != null )
			System.out.println(name);
		for( i=0;i<a.length;i++ )
		{
			a[i].print();
			if( (comma&1) == 1 )
				System.out.printf(",");
		}
		if( (comma&2) == 2 )
			System.out.printf("\n");
	}

	Shift TestShift( int s )
	{
		int j;
		Shift c = this.shiftOne(s);
		Shift b = this.shiftOneByExt(s);
		if( !b.Equ(c) )
		{
			System.out.printf("b = ");
			for( j=0;j<Shift.size;j++ )
				System.out.printf("%08x ",b.buf[j]);
			System.out.printf("\n");
			System.out.printf("c = ");
			for( j=0;j<Shift.size;j++ )
				System.out.printf("%08x ",c.buf[j]);
			System.out.printf("\n");
		}
		
		return b;
	}
	static int TestXk(int k, int s)
	{
		Shift a = new Shift();
		a.setXk(k);
		a = a.TestShift(s);
		return a.isXk();
	}
	private void setBit(int k, int pos)
	{
		
		buf[pos/32] |= (k<<(pos%32));
	}
	static Shift[] msg2shift(byte[] msg)
	{
		int len = msg.length*8/length;
		Shift[] ret = new Shift[len];
		for( int i=0;i<len; i++ )
			ret[i] = new Shift();
		Shift.msg2shift(ret, msg, len*length);
		return ret;
	}
	static int msg2shift( Shift[] info, byte[] msg, int k )
	{
		if( Shift.length*info.length != k )
			return -1;
		int i;
		for( i=0;i<info.length;i++ )
			info[i].clean();
		for( i=0;i<k;i++ )
			info[i/Shift.length].setBit(((int)msg[i/8]>>(i%8))&1,i%Shift.length);
		return Shift.length*info.length;
	}
	public static int shift2msg( Shift[] info, byte[] msg, int k )
	{
		int i;
		for( i=0;i<info.length*Shift.length/8;i++ )
			msg[k+i]=0;
		for( i=0;i<info.length*Shift.length;i++ )
		{
			int pos = i%Shift.length;
			int bit = info[i/Shift.length].buf[pos/32]>>(pos%32);
			bit &= 1;
			msg[k+i/8]|=(bit<<i%8);
		}
		return info.length*Shift.length;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int i,j,k;
		Shift.length = 136;
		Shift.init();
		
		Shift a = new Shift();
		
		for( i=0;i<10000;i++ )
		{
			
			for( j=0;j<Shift.size;j++ )
			{
				a.buf[j] = (int) Math.round(Math.random()*2147483647.);
			}
			a.buf[Shift.size-1] &= Shift.mask;
			
			a.TestShift(i%136);
		
			if( i%1000 == 0 )
				System.out.printf("*");
		}
		System.out.printf("\n");
		
		a.setXk(14);
		k=2;
		for( i=0;i<10000;i++ )
		{
			k = Shift.TestXk(k,i%136);
		}
		
		System.out.printf("\n");
		a.print();
		System.out.printf("\n");
		
		
	}
}
