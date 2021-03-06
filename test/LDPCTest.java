import com.jove.ldpc.*;
import java.io.*;
import java.util.Scanner;

class LDPCTest {

	public static void main( String[] argv ) throws IOException{
		Scanner deg = new Scanner( new File(argv[0]) );
		File ofile = new File(argv[1]);
		double rate = Double.parseDouble(argv[2]);
		int bl = Integer.parseInt(argv[4]);
		int N = Integer.parseInt(argv[3]);
		int poly = Integer.parseInt(argv[5],16);

		int length = deg.nextInt();
		int dg[] = new int[length];
		double df[] = new double[length];
		int i;
		for( i=0;i<length;i++ )
			dg[i]=deg.nextInt();
		for( i=0;i<length;i++ )
			df[i]=deg.nextDouble();

//		int[] dg = {2,3,4,5,8,10,12};
//		double[] df = {0.4717,0.3336,0.0102,0.0426,0.0070,0.0049,0.1300};  
		double[] range = {1.0,1.1};
		Degree Ddeg = new Degree(dg,df);
		LDPCDBItem item = new LDPCDBItem( Ddeg, "", "", range, rate );
		PrintWriter writer = new PrintWriter(new FileWriter(ofile));
		item.logCode(writer,bl,N,poly);
		item.logTask(writer,bl,N,poly);
		QCRS rs = new QCRS(item.getCode());
		for( i=0;i<(rs.N-rs.M)*bl;i+=73 )
			QCRS.testEnc(rs, i);
		rs.reportTask();
		writer.flush();
		writer.close();
	}
}

