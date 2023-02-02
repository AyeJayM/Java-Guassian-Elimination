import java.io.*;
import java.util.*; //Allows for Vector functionality

public class gaussian


{
   
    public static void main(String[] args) throws IOException
    {
        boolean SPP = false;    //This will determine if we use SPP method
        String filename1;       //This will hold the user-provided filename

        //Check if user inputted "-spp"
        if (args[0].equals("-spp") || args[0].equals("--spp"))
        {
            filename1 = args[1];
            SPP = true;
        }
        else
        {
            filename1 = args[0]; //This string holds the inputted filename;
        }


        Scanner scanSys = new Scanner(new File(filename1)); //Set a scanner to file


        
        //To get the first int which states num of entries per row

        int scanned = scanSys.nextInt();


        //This will allow us to determine the number of rows
        BufferedReader reader = new BufferedReader(new FileReader(filename1));
        int totalLines = 0;
        while (reader.readLine() != null)
        {
            totalLines++;
        }
        reader.close();

        double[][] LinSystem = new double[totalLines-2][scanned];   //Create our 2D doubles array
                                                                    // "-1" is because the first int isn't counted
                                                                    //nor the last line which constains constants
        
        double[] con = new double[scanned]; //An array of our constants



        int row = 0;
        int currentLine = 0;

        // Scans doubles from file into a 2D array
        while ( scanSys.hasNextLine() )
        {
            for (int col = 0; col < scanned; col++)
            {
                if (currentLine == totalLines-2)
                {
                    con[col] = scanSys.nextDouble();                 
                }
                else
                {
                    LinSystem[row][col] = scanSys.nextDouble();
                }      
            }
            currentLine++;
            row++;
        }
        scanSys.close();

        
        // This is just to test if it is storing the doubles into an array correctly
        System.out.println("\nNow printing the inputted linear systems...\n");
        for (int PrintRow = 0; PrintRow < row-1; PrintRow++)
        {
            for (int PrintCol = 0; PrintCol < scanned; PrintCol++)
            {
                System.out.print(LinSystem[PrintRow][PrintCol]);
                System.out.print(" ");
            }

            System.out.print("\n");
        }

        System.out.println("\nNow printing the array of constants.\n");
        for (int conCount = 0; conCount < scanned; conCount++)
        {
            System.out.print(con[conCount] + " ");
        }

        System.out.println("\n\nSPP Requested? : " + SPP +"\n");
  

        if (SPP == true)
        {
            SPPGaussian(LinSystem, con, scanned, filename1);
        }
        else 
        {
            NaiveGaussian(LinSystem, con, scanned, filename1);
        }

    }


///////////////////////////////////////////////////////////////////////




    public static void NaiveGaussian(double[][] coeff, double[] cons, int scanned, String filename1)
    {
        double sol[] = new double[scanned];

        FwdElimination(coeff, cons, scanned);
        BackSubst(coeff, cons, sol, scanned, filename1);
    }

    public static void FwdElimination(double[][] coeff, double[] cons, int scanned)
    {

        for (int k = 0; k < scanned-1; k++ )
        {
            for (int i = k + 1; i < scanned; i++)
            {
                double mult = ( coeff[i][k] / coeff[k][k] );

                for (int j = k; j < scanned; j++)
                {
                    coeff[i][j] = coeff[i][j] - mult * coeff[k][j];
                }
                cons[i] = cons[i] - (mult*cons[k]);
            }
        }
    }

    public static void BackSubst(double[][] coeff, double[] cons, double[] sol, int scanned, String filename1)
    {

        //We do "scanned -1" since sol[4] is out of bounds for array of size 4
        sol[scanned-1] = cons[scanned-1]/coeff[scanned-1][scanned-1];

        for (int i = scanned-2; i >= 0; i--)
        {
            double sum = cons[i];

            for (int j = i+1; j < scanned; j++)
            {
                sum = sum - coeff[i][j] * sol[j];
            }
            sol[i] = sum/coeff[i][i];
        }
        
        System.out.println("\nYour solution is...\n");

        for (int result = 0; result < scanned; result++)
        {
            System.out.print(sol[result] + " ");
        }

        System.out.println("\n\nYour solution will be printed to a file of the same name with a \".sol\" extension.\n");

        try {
            SolutionToFile(sol, scanned, filename1);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
    

/////////////////////////////////////////////////////////////////


    public static void SPPGaussian(double[][] coeff, double[] cons, int scanned, String filename1)
    {
        double[] sol = new double[scanned];
        int[] ind = new int[scanned];

        for (int i = 0; i < scanned; i++)
        {
            ind[i] = i;
        }

        SPPFwdElimination(coeff, cons, ind, scanned);
        SPPBackSubst(coeff, cons, sol, ind, scanned, filename1);

    }



    public static void SPPFwdElimination(double[][] coeff, double[] cons, int[] ind, int scanned)
    {
        double[] scaling = new double[scanned]; //Array of scaling factors

        //Initialize index and scaling vectors
        for (int i = 0; i < scanned; i++)
        {
            double smax = 0;

            for(int j = 0; j < scanned-1; j++)
            {
                smax = Math.max( smax , Math.abs(coeff[i][j]) ); //find coefficient with greatest absolute value
            }

            scaling[i] = smax;
        }

        for(int k=0; k < scanned-1 ; k++)
        {
            double rmax = 0;

            int maxInd = k;

            for( int i = k; i < scanned; i++)
            {
                double r = Math.abs(coeff[i][k]/scaling[ind[i]]); //ratio of coefficient to scaling factor

                if (r > rmax)
                {
                    rmax = r;
                    maxInd = i;
                }
            }

            //swap(ind[maxInd], ind[k]);
            //This mimicks that functionality.

            int temp = ind[maxInd];

            ind[maxInd] = ind[k];

            ind[k] = temp;


            for(int i = k+1 ;i<scanned; i++)
            {
                double mult = coeff[ind[i]][k] / coeff[ind[k]][k];

                for(int j = k; j<scanned; j++)
                {
                    coeff[ind[i]][j] = coeff[ind[i]][j] - mult*coeff[ind[k]][j];
                }

                cons[ind[i]] = cons[ind[i]] - mult*cons[ind[k]];

            }
            
        }

    }




    public static void SPPBackSubst(double[][] coeff, double[] cons, double[] sol, int[] ind, int scanned, String filename1)
    {
        sol[scanned-1] = cons[ind[scanned-1]] / coeff[ind[scanned-1]][scanned-1];

        for(int i = scanned - 2; i >= 0; i--)
        {
            double sum = cons[ind[i]];

            for(int j=i+1; j<scanned; j++)
            {
                sum = sum - coeff[ind[i]][j] * sol[j];
            }

            sol[i] = sum/coeff[ind[i]][i];
        }

        System.out.println("\nYour solution is...\n");

        for(int printResults = 0; printResults < scanned; printResults++)
        {
            System.out.print(sol[printResults] + " ");
        }

        System.out.println("\n");
        System.out.println("\nYour solution will be printed to a file of the same name with a \".sol\" extension.\n");
        
        try {
            SolutionToFile(sol, scanned, filename1);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /////////////////////////////////////////////////////////////


    public static void SolutionToFile(double[] solution, int scanned, String filename1) throws IOException
    {
        //We will first need to remove the ".lin" extension from the filename
        String filename2 = ( filename1.substring(0, filename1.length() - 3) ) + "sol";

        BufferedWriter writer;
        writer = new BufferedWriter(new FileWriter(filename2));

        for (int i = 0; i < scanned; i++)
        {
            writer.write(String.valueOf(solution[i]));
            writer.write(" ");
        }

        writer.close();;
    }

}

