package com.jove.ldpc;

import java.util.LinkedList;
import java.util.List;

public class BlockDec {

	class stCodePair
	{
		int pos;
		int shift;
	};
	class stCodeBlock
	{
		int row;
		int col;
		int shift;
		float[] pr0d1;
		float[] pr0d1B;
		void backup()
		{
			for(int i=0;i<pr0d1.length;i++ )
				pr0d1B[i] = pr0d1[i];
		}
	};
	class stCodeCol
	{
		int cw;
		stCodeBlock[] blk;
		void backup()
		{
			for( int i=0;i<cw;i++ )
				blk[i].backup();
		}
	};
	int col;
	int row;
	int bl;
	int len;
	int[][] codeTable;
	int intp;

	public BlockDec( int[][] codeTable, int blockLength )
	{
		bl = blockLength;
		col = codeTable[0].length;
		row = codeTable.length;
		this.codeTable = codeTable;
		buildCode(codeTable);
	}

	public int dec( float[] in, float std_de, int cnt, byte[] out )
	{
		int i,j,k;
		float[][] tempp;
		int[] check=new int[1];
		initDec( in, std_de );

		for( i=0;i<col;i++ )
			setToOne( oldprp[i] );
		for( i=0;i<row;i++ )
		{
			for( j=0;j<CodeCol[i].cw;j++ )
				setToOne( CodeCol[i].blk[j].pr0d1 );
		}	
		for( i=0;i<cnt;i++ )
		{
			for( j=0;j<row;j++ )
				CodeCol[j].backup();
			for( j=0;j<col;j++ )
				setToOne( newprp[j] );
			check[0]=0;
			decOnce(check);
			tempp=oldprp;
			oldprp=newprp;
			newprp=tempp;
			if(check[0]==0)
				break;
		}
		for(j=0;j<out.length;j++ )
			out[j]=0;
		for( j=0;j<col;j++ )
			for( k=0;k<bl;k++ )
			{
				int r=(oldprp[j][k]*prz[j][k]>1f)?0:1;
				int pos = j*bl+k;
				out[pos/8]|=(r<<(pos%8));
			}
		return 0;
	}

	float[] work;
	float[] dlp;
	float[] tbuf;
	float[][] newprp;
	float[][] oldprp;
	float[][] prz;
	stCodePair[] CodePair;
	stCodeCol[] CodeCol;
	class stCodeRow
	{
		List<stCodeBlock> list;
		stCodeRow()
		{
			list = new LinkedList<stCodeBlock>();
		}
	}
	stCodeRow[] CodeRow;
	
	int initDec( float[] in, float std_de )
	{
		int i;
		for( i=0;i<col;i++ )
		{
			float[] p=prz[i];
			for(int j=0;j<bl;j++ )
				p[j]=vflimit((float)Math.exp(2.f*(float)in[i*bl+j]/((float)std_de*(float)std_de)));	
		}
		return 0;

	}

	int weightCode( int[][] code )
	{
		int i;
		int cw,decp;
		decp=0;
		for( i=0;i<row;i++ )
		{
			cw=buildCodePair( code[i] );
			decp+=cw;
		}
		return decp;
	}
	float[] allocBL()
	{
		return new float[bl];
	}
	float[][] allocWork()
	{
		float[][] ret = new float[col][];
		for( int i=0; i<col; i++ )
			ret[i]=allocBL();
		return ret;
	}
	int buildCode( int[][] code )
	{
		int i,j;
		int cw;
		CodePair=new stCodePair[col];
		CodeCol=new stCodeCol[row];
		CodeRow=new stCodeRow[col];
		work=allocBL();	
		dlp=allocBL();	
		tbuf=allocBL();	
		prz= allocWork();
		oldprp = allocWork();
		newprp = allocWork();
		for( i=0;i<col;i++ )
			CodeRow[i] = new stCodeRow();
		for( i=0;i<row;i++ )
		{
			cw=buildCodePair( code[i] );
			CodeCol[i] = new stCodeCol();
			CodeCol[i].blk=new stCodeBlock[cw];
			CodeCol[i].cw=cw;
			for( j=0;j<cw;j++ )
			{
				CodeCol[i].blk[j] = new stCodeBlock();
				CodeCol[i].blk[j].row=i;
				CodeCol[i].blk[j].col=CodePair[j].pos;
				CodeCol[i].blk[j].shift=CodePair[j].shift;
				CodeCol[i].blk[j].pr0d1=new float[bl];
				CodeCol[i].blk[j].pr0d1B=new float[bl];
				CodeRow[CodePair[j].pos].list.add(CodeCol[i].blk[j]);
			}
		}
		return 0;
	}
	int buildCodePair( int[] c )
	{
		int cw=0;
		int i;
		for( i=0;i<col;i++ )
		{
			if( c[i]!=bl )
			{
				CodePair[cw] = new stCodePair();
				CodePair[cw].pos=i;
				CodePair[cw].shift=c[i];
				if( c[i]<0 )
					System.out.println("Error Shift "+c[i]);
				cw++;		
			}
		}
		return cw;
	}
	float vfmax( float a, float b )
	{
		return (a>b)? a:b;
	}
	float vfmin( float a, float b )
	{
		return (a<b)? a:b;
	}
	static float vfESPL = (float)1e-20;
	static float vfESPH = (float)1e20;
	float vflimit( float a )
	{
		return vfmin(vfmax(a,vfESPL),vfESPH); 
	}
	int blockMove(float[] dec, float[] src, int shift)
	{
		int i;
		for(i=0;i<bl;i++)
		{
			dec[i]=src[(i+shift)%bl];
		}
		return bl;
	}

	int decBlockKernel( int pos, int[] check )
	{
		float v0,v1;
		float[] pa,pb;
		int num,i,j;
		int cw = CodeCol[pos].cw;
		stCodeBlock[] blk = CodeCol[pos].blk;
		float[][] pr1m0 = new float[cw][];
		
		setToOne(dlp);

		for( i=0;i<cw;i++ ) // BUG dlp only need one per Col
		{
			pr1m0[i] = new float[bl];
			
			blockMove(tbuf, oldprp[blk[i].col], blk[i].shift );
			blockMove(work, prz[blk[i].col], blk[i].shift );

			for( num=0;num<bl;num++ )
			{	
				v0=work[num]*tbuf[num]/blk[i].pr0d1[num];		
				v1=1f-2f/(v0+1f);	/*LR=(PR-1)/(PR+1)? (1-PR)/(1+PR)*/
				if(v1<0f)
				{
//					System.out.println(num+" "+work[num]+" "+tbuf[num]+" "+blk[i].pr0d1[num]);
//					float v2 = 1f;
//					for( stCodeBlock s : CodeRow[blk[i].col].list )
//					{
//						int p = (num - s.shift+blk[i].shift+bl)%bl;
//						v2*=s.pr0d1B[p];
//						System.out.println("p0/p1 blk["+s.col+","+s.row+","+s.shift+"].pr0d1["+p+"]="+s.pr0d1B[p]);
//					}
//					if( v2!=tbuf[num] )
//						System.out.println("Error: "+tbuf[num]+" "+v2);
				}
				dlp[num]=dlp[num]*v1;
				pr1m0[i][num]=v1;				
			}
		}
		for( num=0;num<bl;num++ )
		{
			if( dlp[num]<0 )
			{
				check[0]=1;
//				for( i=0;i<cw;i++ )
//					System.out.println("p1-p0 blk["+blk[i].col+"].pr1m0["+num+"]="+pr1m0[i][num]);
			}
		}
		for( i=0;i<cw;i++ )
		{
			for( num=0;num<bl;num++ )
			{	
				v0=dlp[num]/pr1m0[i][num];
				v1=(1f+v0)/(1f-v0);			/*PR=(1+LR)/(1-LR) (1-LR)/(1+LR) */
				blk[i].pr0d1[num]=vflimit(v1);	
			}
			blockMove( work,blk[i].pr0d1,(bl-blk[i].shift+bl)%bl );
			pa=work;
			pb=newprp[blk[i].col];
			for( num=0;num<bl;num++ )
			{
				pb[num]*=pa[num];
			}		
		}
		
		return 0;
	}
	int decOnce(int[] check)
	{
		int i;
		for( i=0;i<row;i++ )
		{
			decBlockKernel(i,check);
		}
		return 0;
	}

	void setToOne( float[] buf )
	{
		int i;
		for( i=0;i<bl;i++ )
			buf[i]=1f;
	}



	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int[][] rs = {
				{  228,  234,  256,  227,  253,  256,  256,  248,  256,  256,  256,    0,  256,  256,  256,  256,  256,  256,  256,  256, },
				{  256,   17,   13,  256,  247,  256,   15,  256,  256,  249,  256,  256,    0,  256,  256,  256,  256,  256,  256,  256, },
				{  256,  227,  252,  241,  256,  256,  256,  230,    2,  256,  256,  256,  256,    0,  256,  256,  256,  256,  256,  256, },
				{  256,   24,   29,  256,  256,    7,  256,   29,  256,  256,   30,  256,   20,  256,    0,  256,  256,  256,  256,  256, },
				{  256,  256,  245,  231,  256,  256,  256,  256,  232,  256,    0,  251,  256,  256,  256,    0,  256,  256,  256,  256, },
				{  256,  252,  256,    8,  256,   15,  251,  256,  256,  256,    7,  256,  256,   12,  256,  256,    0,  256,  256,  256, },
				{  256,  254,  256,   11,  256,    8,  256,  256,  256,   11,  240,  256,  256,  256,  256,  256,  256,    0,  256,  256, },
				{  236,  256,  234,    1,  256,  256,  245,  256,  256,  256,  252,  256,  256,  256,  256,  256,  256,  256,    0,  256, },
				{  256,  256,   15,  244,  256,  256,  256,  256,  256,  256,    4,  256,  256,  256,  256,  250,    7,   19,  256,    0, },
				{  256,  250,  235,  256,    3,  256,  256,  256,  252,  256,    1,  256,  256,  256,    0,  256,  256,  256,  240,  239, },
				};
		QCRS enc = new QCRS(rs);
		enc.encoderInit();
		float sigma = 0.80f;
		byte[] msg = new byte[rs[0].length*256/8];
		byte[] out = new byte[rs[0].length*256/8];
		int i,j;
		Shift[] so;
		int infoLen = rs[0].length - rs.length;
		infoLen *= 256;
		for( i=0;i<infoLen/8;i++ )
			msg[i] = (byte)(i&0xff);
		msg[31] = (byte)(0xaa&0xff);
		
		int ret = Shift.msg2shift(enc.info, msg, infoLen);
		so = enc.enc();
		enc.verify(15);
		Shift.shift2msg(so, msg, 0 );
		AWGN channel = new AWGN(sigma);
		float[] recv = channel.channel(msg);
		float n = AWGN.calcN(recv);
		System.out.println("Noise :"+n);
		BlockDec dec = new BlockDec(rs,256);
		dec.dec(recv, sigma, 100, out);
		so = Shift.msg2shift(out);
		Shift.printArray(so, "Result ", 3);
		for( i=0;i<out.length;i++ )
			if( out[i]!=msg[i] )
				System.out.println(i+"  Error "+(int)out[i]+" "+(int)msg[i]);
	}

}
