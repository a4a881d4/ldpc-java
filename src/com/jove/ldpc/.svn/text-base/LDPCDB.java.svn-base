package com.jove.ldpc;

import java.util.LinkedList;
import java.util.List;

public class LDPCDB {
		
	public List<LDPCDBItem> db;

	void append( LDPCDBItem aCode )
	{
		db.add(aCode);
	}
	
	void append( int dg[], double df[], String comm, String url, double[] Eb, double rate )
	{
		Degree deg = new Degree(dg,df);
		LDPCDBItem item = new LDPCDBItem( deg, comm, url, Eb, rate );
		append(item);
	}
	
	static double[] range(double from, double to, double step )
	{
		int len = (int)Math.floor((to-from)/step)+1;
		double[] ret= new double[len];
		for( int i=0;i<len;i++ )
			ret[i] = i*step+from;
		return ret;
	}
	public LDPCDB()
	{
		db = new LinkedList<LDPCDBItem>();
		append( new int[]{2,3,4,5,8,10,12},
				new double[]{ 0.4717,0.3336,0.0102,0.0426,0.0070,0.0049,0.1300 },
				"Design of Capacity-Approaching Irregular Low-Density Parity-Check Codes sigma = 0.9580",
				"",
				LDPCDB.range(1.5, 0.8, -0.1),
				0.5
				);
		append( new int[]{2,3,12},
				new double[]{ 0.750000,0.166667,0.083333333 },
				"Unknow rate 0.125",
				"",
				LDPCDB.range(1.5, 0.8, -0.1),
				0.125
				);
		append( new int[]{2,3,5},
				new double[]{ 0.521814,0.271293,0.206893 },
				"Unknow rate 0.5",
				"",
				LDPCDB.range(1.5, 0.8, -0.1),
				0.5
				);
		append( new int[] {2,3,4,11},
		    new double[] {0.450013,0.370771,0.0307238,0.1484922},
		    "rate 0.5",
		    "",
		    LDPCDB.range(1.5, 0.8, -0.1),
        0.5
		    );
		append( new int[] {2,3,4,5,8,10,12},
        new double[] {0.4717,0.3336,0.0102,0.0426,0.0070,0.0049,0.1300},
        "rate 0.25",
        "",
        LDPCDB.range(1.5, 0.8, -0.1),
        0.25
        );
    
//		append( new int[]{2,3,4,5,6,7,8,9,10,11,12},
//				new double[]{0.242187,0.189400,0.061625,0.042975,0.063281,0.064757,0.004461,0.148330,0.007488,0.076095,0.099402},
//				"7,8,0.308409,0.691591 0.930000 < sigma < 0.935000",
//				"http://sonic.newcastle.edu.au/ldpc/lopt/results/a738431dc",
//				0.5
//				);
//		append( new int[]{2,3,4,12},
//				new double[]{0.215176,0.180078,0.384898,0.219849},
//				"7,8,0.900454,0.099546 0.925000 < sigma < 0.930000",
//				"http://sonic.newcastle.edu.au/ldpc/lopt/results/ae585316b",
//				0.5
//				);
	}
	public static void main(String[] args) 
	{
		LDPCDB db = new LDPCDB();
		LDPCDBItem x = db.db.get(3);
		int bl = 96;
		int N = 48;
		int poly = 0x83;
		int SIM = 0;
		int i;
		if( x.getCode(N, bl, poly, null) != null)
		{
			QCRS enc = new QCRS(x.code);
			enc.encoderInit();
			byte[] msg = new byte[N*bl/8];
			byte[] out = new byte[N*bl/8];
			for( i=0;i<(enc.N-enc.M)*enc.bl/8;i++ )
				msg[i] = (byte)Math.round(Math.random()*255.);
//			for( i=0;i<(enc.N-enc.M)*enc.bl/8;i++ )
//				msg[i] = (byte)(0&0xff);
//			msg[0]=2;
			enc.setInfo(msg);
			Shift[] so = enc.enc();
			enc.verify(15);
			Shift.printArray(so, " enc ", 3);
			Shift.shift2msg(so, msg, 0 );
			BlockDec dec = new BlockDec(x.code,bl);
			float[] sigma = x.sigma();
			for( int j=0;j<sigma.length;j++  )
			{
				AWGN channel = new AWGN(sigma[j]);
				int err=0;
				for( int cnt=0;cnt<SIM;cnt++ )
				{
					float[] recv = channel.channel(msg);
					dec.dec(recv, sigma[j], 100, out);
					//Shift[] s1 = Shift.msg2shift(out);
					//Shift.printArray(s1, " dec ", 3);
					for( i=0;i<out.length;i++ )
						if( out[i]!=msg[i] )
						{
							System.out.println(x.EbN0Range[j] + "  "+ i+"  Error "+(int)out[i]+" "+(int)msg[i]);
							err++;
							break;
						}
				}
				System.out.println("EbN0 "+x.EbN0Range[j]+"Blk Err Rate "+(double)err/(double)SIM);
			}	
		}
		int[][] ttu = x.code;
		System.out.printf("%d %d %d\n",N,ttu.length,bl);
    for( i=0;i<ttu.length;i++ )
    {
      for( int j=0;j<N;j++ )
        System.out.printf("%4d ",ttu[i][j]);
      System.out.printf("\n");
    }
	}
}
