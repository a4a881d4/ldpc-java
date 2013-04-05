import com.jove.ldpc.*;
import java.io.*;

class LDPCTest {

	public static void main( String[] argv ) throws IOException{
		File file = new File(argv[0]);
		int[] dg = {2,3,4,5,8,10,12};
		double[] df = {0.4717,0.3336,0.0102,0.0426,0.0070,0.0049,0.1300};  
		double[] range = {1.0,1.1};
		Degree deg = new Degree(dg,df);
		LDPCDBItem item = new LDPCDBItem( deg, "", "", range, 0.25 );
		PrintWriter writer = new PrintWriter(new FileWriter(file));
		item.logCode(writer,72,28,0x83);
		writer.flush();
		writer.close();
	}
}

